import configparser

# This is used by the shasta scripts to access the config file.

# Read the config file.
def getConfig():
    config = configparser.ConfigParser()
    configFileName = 'shasta.conf'
    if not config.read(configFileName):
        print('Error reading config file %s.' % configFileName)
    return config


