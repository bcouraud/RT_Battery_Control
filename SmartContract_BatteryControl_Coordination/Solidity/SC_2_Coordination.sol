// This Smart Contract opens a negotiation place (Market Place): First, it will collect the bids from buyers and sellers. It will store the bids.
// Then, it will proceed to the paymen of the selected bids, where the selection is done outside the Smart Contract, by an external agent (DSO)

pragma solidity >=0.4.21 <0.7.0;

contract BatteryCoordination {
    uint8 private householdsCount;
    uint private operationalDeposit; // an amount of money used to pay all operation fees
    uint private aggregatedSoC;
    uint private energycommitted;
    uint private energyalreadyexported;
//    uint private SoC; // an amount of money used to pay all operation fees
    // The market place stores (with a hash) the bid (price and quantity, we also need the weights (or Beta) and other informations to be added), and the type of the agent (Buyer (1) or Seller (0))
 //   mapping (address => uint) private bidPrice;
 //   mapping (address => uint) private bidQuantity;
 //   mapping (address => uint) private bidWeight;
 //   mapping (address => uint) public agentType;
 //   mapping (address => uint) private depositAmount;
    mapping (address => uint) private depositAmount;
    mapping (address => uint) private SoC;
    mapping (address => uint) private EnergyExported;
    mapping (address => uint) private WeightImport;
    mapping (address => uint) private WeightExport;
    //If we use state machine to coordinate all steps
    enum State { DEPOSIT_OPEN, DATA_COLLECTION, COORDINATION, COMPLETE, FINISHED}
    State public currState;

    // The market place is owned by the DSO Agent:
    address payable public AggregatorAgent;
    // We also need to add a time: Times are either absolute unix timestamps (seconds since 1970-01-01) or time periods in seconds.
    uint public EndTime; // not used yet

       // We need  a boolean to state the stage of the coordination phase (open or not for receiving new offers/bids) Set to true at the end, disallows any change. By default initialized to `false`.
    bool ended;

    //  We create a modifier that could be used later to restrict some functions to only the operator. not needed for simple use
    modifier onlyOperator() {
        require(
            msg.sender == AggregatorAgent,
            "Only Aggregator can call this."
        );
        _;
    }

    // Constructor function to initialize the MarketPlace. It is "payable" so it can receive initial funding to cover up some mispayment
    constructor() public payable {
       // require(msg.value == 10 ether, "10 ether initial funding required");
        /* Set the owner to the creator of this contract */
        AggregatorAgent = msg.sender;
        // initialization of variables
        householdsCount = 0;
        aggregatedSoC = 0;
    }


    // used to refill the amount of the operational deposit
    function operationaldeposit() public payable onlyOperator returns (uint) {
        operationalDeposit = operationalDeposit+msg.value;
         return 1;
    }

    function getaggregatedSoc() public payable onlyOperator returns (uint) {
         return aggregatedSoC;
    }
     function resettargets() public payable onlyOperator returns (uint) {
        energycommitted=100;
        energyalreadyexported=0;

         return 1;
    }   

    ///  Agents submit their information
    function submitHouseholdData(uint _soc, uint _energyexported, uint _agenttype) public payable  returns (uint) {
         // if we use state machine, require(currState == State.DATA_COLLECTION, "Cannot confirm data submission");
       if (_agenttype == 1) {
            householdsCount++;}
        else  {

        }
        // We update all values of this agent
        aggregatedSoC = aggregatedSoC - SoC[msg.sender] + _soc;
        energyalreadyexported = energyalreadyexported - EnergyExported[msg.sender]+ _energyexported;
        SoC[msg.sender] = _soc;
        EnergyExported[msg.sender] = _energyexported;
        if (_energyexported >= 100) {
            WeightExport[msg.sender] = _soc*1000/aggregatedSoC;}
        else  {
            WeightExport[msg.sender] = 0;
        }
        WeightImport[msg.sender] = 0;
        depositAmount[msg.sender] = msg.value;
        return WeightExport[msg.sender];
    }

/*  function getbalance(address agentaddress) public view returns(uint) {
        uint out;
        out = depositAmount[agentaddress];
        return out;
    }
  function getrealbalance(address agentaddress) public view returns(uint) {
        uint out;
            if (agentaddress!=AggregatorAgent){
                out = agentaddress.balance;
            }
            else
            {
                out = AggregatorAgent.balance;
            }
        return out;
    }*/

    function ExtractWeights(address agentaddress) public view returns(uint[2] memory) {  // we store the
    // return into memory, not storage, as we do not ned it outside of the function
        require(msg.sender == AggregatorAgent,"Only DSO/Operator can call this.");
        // if we use state machine, require(currState == State.BID_DEPOSIT_OVER, "Cannot confirm Bid deposit");
        uint[2] memory array;

            array = [WeightImport[agentaddress], WeightExport[agentaddress]];
            return array;
         //currState = State.COORDINATION;

    }

// Once everything is finished, the DSO/operator close the negotiation by redistributing the money that is in the accounts
  function close(address payable agentaddress) external payable onlyOperator {
        // if we use state machine, require(currState == State.COMPLETE, "Cannot confirm delivery");
            // selleraddress.transfer(depositAmount[selleraddress]); // we transfer the money to the seller
            if (agentaddress!=AggregatorAgent){
            agentaddress.transfer(depositAmount[agentaddress]); // we transfer the money to the agent
            }
            else
            {
                AggregatorAgent.transfer(operationalDeposit);
            }
        //currState = State.COMPLETE;
    }

    // Once everything is finished, the aggregator closes the steps by redistributing the money to the accounts
  function closecontract() external payable onlyOperator {
            uint fee;
            fee = 0;  // we define is a fee of operational cost that should be removed. 0 for private blockchain.
                msg.sender.transfer(operationalDeposit-fee);
        //currState = State.FINISHED;
    }


}

