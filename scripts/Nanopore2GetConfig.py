import configparser

# This is used by the Nanopore2 scripts to access the config file.

# Read the config file.
def getConfig():
    config = configparser.ConfigParser()
    configFileName = 'Nanopore2.conf'
    if not config.read(configFileName):
        print('Error reading config file %s.' % configFileName)
    return config


