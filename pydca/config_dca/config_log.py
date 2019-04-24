
"""A sample configuration for logging. This is a suggestion and you are free to
chose your own preferred way of configuring loggers.
"""

LOGGING_CONFIG = {
    'version':1,
    'disable_existing_loggers':False,
    'formatters':{
        'verbose':{
            'format':'%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s'
        },
        'simple':{
            'format':'%(levelname)s %(message)s'
        },
    },
    'handlers':{
        'console':{
		    'level':'INFO',
			'class':'logging.StreamHandler',
            'formatter':'verbose',
		},
	},
	'loggers':{
		'':{
			'handlers':['console'],
			'level':'DEBUG',
			'propagate':True,
		},
	},

}

class ConsoleColor:
    """Defines colors for logging messages.

    Attributes
    ----------
        nocolor :  str
            Disable coloring of terminal while logging messages.
        green : str
            This color is used for logging level INFO.
        yellow : str
            This color is used for logging level WARNING.
        red : str
            This color is used for logging level ERROR.
    """
    nocolor = '\033[0;0m'
    green = '\033[0;32m'
    yellow = '\033[0;33m'
    red = '\033[0;31m'
