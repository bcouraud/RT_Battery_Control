# """ Before running this, you need to install Ganache to run a local Blockchain. Otherwise, use infura to access a node. """
# Here, addresses match the local Ganache Blockchain.
# Then, You should (maybe) install solc (open a powershell terminal or a cmd, and run "pip install solc")
# Then, You should (definitely) install solcx (open a powershell terminal or a cmd, and run "pip install py-solc-x")  https://pypi.org/project/py-solc-x/
# Finally, you should change the paths that are listed below to locate the solidity SC (here, Greeting.sol).
# import json
import time
import pprint
from web3 import Web3
from solcx import compile_source, compile_files
from solcx import install_solc
install_solc('v0.5.3')
import numpy as np
import pandas as pd

#########################  To Run the interaction with Matlab ############################
# First, follow the instructions below: https://uk.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html 
# 1. Make sure you have Python in your PATH.
# 2. Find the Matlab root folder. You can use the matlabroot command within Matlab to find it.
# 3. In the windows Command line ("cmd" opened with admin rights) Go to the Matlab root folder in the command line by typing cd "matlabroot\extern\engines\python" (In Windows)
# 4. Type in python setup.py install 
# 6. In matlab (opened with admin rights ) type in :cd (fullfile(matlabroot,'extern','engines','pytpythohon'))
#                                                      system('python setup.py install')
# Then, uncomment the 2 following lines:
""" import matlab.engine
matlab_eng = matlab.engine.start_matlab() """




def compile_source_file(file_path):
   with open(file_path, 'r') as f:
      source = f.read()

   return compile_source(source)


def deploy_contract(w3, contract_interface):
    tx_hash = w3.eth.contract(
        abi=contract_interface['abi'],
        bytecode=contract_interface['bin']).constructor().transact()
#functions.transact() executes the specified function by sending a new public transaction.
    address = w3.eth.getTransactionReceipt(tx_hash)['contractAddress']
    return address

ether = 10**18 # 1 ether = 1000000000000000000 wei


# Connection to he Local Ganache Blockchain
ganache_url = "HTTP://127.0.0.1:7545"
w3 = Web3(Web3.HTTPProvider(ganache_url))
print(w3.isConnected())
print(w3.eth.blockNumber)

# We define the agent's accounts
# The aggregator:
w3.eth.defaultAccount = w3.eth.accounts[0]
AggregatorAccount = w3.eth.defaultAccount
# The Agent 1:
Household1Account = w3.eth.accounts[1]
# Agent 2:
Household2Account = w3.eth.accounts[2]

# Compile the contract
contract_source_path = 'c:/Users/bc111/OneDrive - Heriot-Watt University/Smart-Grids/Blockchain/Smart_Contracts/SmartContractVSCodeBatteryControl/Solidity/SC_2_Coordination.sol'
compiled_sol = compile_source_file('c:/Users/bc111/OneDrive - Heriot-Watt University/Smart-Grids/Blockchain/Smart_Contracts/SmartContractVSCodeBatteryControl/Solidity/SC_2_Coordination.sol')
contract_id, contract_interface = compiled_sol.popitem()

# retrieve the compilation results (abi and bytecode)
abi = contract_interface['abi']
bytecode = contract_interface['bin']
# print(abi)
# Deployment of the contract
address = deploy_contract(w3, contract_interface)
contract = w3.eth.contract(
    address = address,
    abi = abi
)

print("Deployed {0} to: {1}\n".format(contract_id, address))

tx = {
      'from': AggregatorAccount,
      'to': address,
      'value': w3.toWei(0,'ether'),
      'gas': 1000000,
      'gasPrice': w3.toWei('5','gwei'),
}
tx_hash = contract.functions.operationaldeposit().transact(tx)
receipt = w3.eth.waitForTransactionReceipt(tx_hash)
#print ("Gas Price: {0}\n".format(contract.getGasPrice()))
print(w3.eth.gasPrice)

gas_estimate = contract.functions.submitHouseholdData(10, 10, 1).estimateGas() #submitHouseholdData the function defined in the SC_1_Bid_and_Payment.sol smart contract - 
# arg(1) is the SoC
# arg(2) is the energy export
# arg(3) is the agent type 
print("Gas estimate to transact with submitHouseholdData: {0}\n".format(gas_estimate))

if gas_estimate < 200000:
    # The seller submits his offer:
  print("Submitting Data to Contract\n")
  tx_hash = contract.functions.submitHouseholdData(550, 10, 1).transact({'from': Household1Account}) # the 
  print("Aggregated SoC : {0}\n".format(contract.functions.getaggregatedSoc().call()))

  #functions.transact() executes the specified function by sending a new public transaction.
  receipt = w3.eth.waitForTransactionReceipt(tx_hash)

  tx = {
        'from': Household2Account,
        'to': address,
        'value': w3.toWei(0,'ether'),
        'gas': 1000000,
        'gasPrice': w3.toWei('5','gwei'),
}
  tx_hash = contract.functions.submitHouseholdData(10000, 15, 1).transact(tx)
  receipt = w3.eth.waitForTransactionReceipt(tx_hash)
  print(receipt)
  print("Aggregated SoC : {0}\n".format(contract.functions.getaggregatedSoc().call()))

else:
  print("Gas cost exceeds 200000")

print("ok")
tx_hash = contract.functions.submitHouseholdData(550, 10, 1).transact({'from': Household1Account}) # the 
receipt = w3.eth.waitForTransactionReceipt(tx_hash)
# we check that the aggregator retrieve the agents's data and computes the weights
print("Household2's Weights [WeightImport*1000, WeightExport*1000] : {0}\n".format(contract.functions.ExtractWeights(Household2Account).call()))
print("Household1's Weights [WeightImport*1000, WeightExport*1000] : {0}\n".format(contract.functions.ExtractWeights(Household1Account).call()))

SoC= pd.read_csv("SoC.csv")
Energy= pd.read_csv("Energy.csv")

# iteration over all agents
Numberofagents = 70
for i in range(Numberofagents):
  gas_estimate = contract.functions.submitHouseholdData(SoC[i], Energy[i], 1).estimateGas() #submitHouseholdData the function defined in the SC_1_Bid_and_Payment.sol smart contract - 
  print("Gas estimate to transact with submitHouseholdData: {0}\n".format(gas_estimate))

  if gas_estimate < 200000:
    print("Submitting Data to Contract\n")
    tx_hash = contract.functions.submitHouseholdData(SoC[i], Energy[i], 1).transact({'from': w3.eth.accounts[i]}) # the 
    print("Aggregated SoC : {0}\n".format(contract.functions.getaggregatedSoc().call()))

    #functions.transact() executes the specified function by sending a new public transaction.
    receipt = w3.eth.waitForTransactionReceipt(tx_hash)

    tx = {
          'from': w3.eth.accounts[i],
          'to': address,
          'value': w3.toWei(0,'ether'),
          'gas': 1000000,
          'gasPrice': w3.toWei('5','gwei'),
  }
    receipt = w3.eth.waitForTransactionReceipt(tx_hash)
    print(receipt)
    print("Aggregated SoC : {0}\n".format(contract.functions.getaggregatedSoc().call()))
    print("Weights : {0}\n".format(contract.functions.ExtractWeights(w3.eth.accounts[i]).call()))

  else:
    print("Gas cost exceeds 200000")

  
print("over")


