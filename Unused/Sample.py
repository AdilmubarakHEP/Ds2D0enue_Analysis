import basf2 as b2
import modularAnalysis as ma

my_path = b2.create_path()

ma.inputMdst(environmentType='default',
             filename="/home/belle2/amubarak/C00-Generation/ccbarDs+EventGeneration.root",
             path=my_path)

ma.fillParticleList("e+:Reco", cut='electronID > 0.5 and abs(dr) < 1 and abs(dz) < 3', path=my_path)

ma.fillConvertedPhotonsList(decayString='gamma:V0_array -> e+ e-',
                            cut='', 
                            path=my_path)
ma.fillParticleList("e+:V0_array",'isDescendantOfList(gamma:V0_array,1)==1', path=my_path)

output_file = 'Reconstruction.root'

ma.variablesToNtuple('e+:Reco', ['sourceObjectIsInList(e+:V0_array)','isInList(e+:V0_array)'],
                     filename=output_file, treename='etree', path=my_path) 

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)