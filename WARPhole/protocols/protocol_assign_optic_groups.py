"""
Create an xml_dir
Launch the script to split into optic groups every x micrographs
The script generates a star files
Parse the star file
Keep a dictionary of micrographs and optics groups that are parsed from the star file
update the optic groups number (add depending on the number of iterations)
Use the star file to add the optic groups to the open particle set (create a new copy)
Re-use the code from:
https://scipion-em.github.io/docs/_modules/relion/protocols/protocol_assign_optic_groups.html#ProtRelionAssignOpticsGroup
"""

import EPU_Group_AFIS #Placeholder, check the correct import
import pyworkflow.utils as pwutils
import os
from datetime import datetime as dt
import emtable
from pyworkflow.object import Integer
from relion.convert.convert31 import OpticsGroups, getPixelSizeLabel

def assignEPUGroupAFIS(micSet,XMLpath):
    importXmlFiles(micSet,XMLpath)
    runAFISscript(XMLpath,starFile)
    micDict = readOpticsGroupStarFile(starFile)
    return(micDict)

def importXmlFiles(micSet,XMLpath):
    micPaths = [mic.filename() for mic in micSet]
    xmlPaths = [pwutils.path.replaceBaseExt(p, 'xml') for p in micPaths]
    subfolder = dt.now().strftime("%H%M%S")
    pwutils.path.makeFilePath(os.path.join(subfolder,path) for path in xmlPaths])
    for p in xmlPaths:
        pwutils.path.copyFile(os.path.join(XMLpath,p), os.path.join(subfolder,p))

def runAFISscript(XMLpath,outputStarFile):
    pass ## TODO: runAFISscript

def readOpticsGroupStarFile(starFile):
    og = OpticsGroups.fromStar(starFile)
    micTable = emtable.Table(fileName=inputStar,tableName='micrographs')
    micDict = {row.rlnMicrographName: row.rlnOpticsGroup for row in micTable}
    return(micDict)

def shiftMicDict(micDict,i):
    for key, value in micDict.iteritems():
        value = value + i
    return(micDict)

def addOpticsGroup(micSet,micDict):
    for mic in micSet:
        ogNumber = micDict.get(mic.filename(),1)
        if not hasattr(mic, '_rlnOpticsGroup'):
            mic._rlnOpticsGroup = Integer()
        mic._rlnOpticsGroup.set(Integer(ogNumber))

############################################
# In protocol_assign_optic_groups
############################################
micDict = assignEPUGroupAFIS(newMicSet,XMLpath)
micDict = shiftMicDict(micDict,i)
self.micDict = {**self.micDict, **micDict}
addOpticsGroup(micSet,self.micDict)

'''
EPU_Group_AFIS.py [-h] --xml_dir XML_DIR [--n_clusters N_CLUSTERS]
                         [--apix APIX] [--mtf_fn MTF_FN] [--voltage VOLTAGE]
                         [--cs CS] [--q0 Q0] [--ftype FTYPE]
                         [--movie_dir MOVIE_DIR] [--output_fn OUTPUT_FN]
                         [--quiet]
'''
