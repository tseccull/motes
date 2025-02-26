#!/usr/bin/env python3


def motes_logs(cwd, verbose):
    '''
    Define the logging configuration with a dict.
    
    Args:
     -- cwd (str)
          Path to the current working directory.
     -- verbose (bool)
          If true, switches on verbose mode for the logger, where all
          logs are also printed to the terminal as well as send to a log
          file.
          
    Returns:
     -- logging_schema (dict)
          Dictionary containing the logging configuration for the motes
          logger.
    '''
    
    if verbose:
        console_level = "INFO"
    else:
        console_level = "WARNING"
	
	# Logging formatter
    motes_formatter = {
        "format" : "%(asctime)s | %(levelname)8s | %(module)16s | %(funcName)18s | %(message)s",
        "datefmt"   : "%Y-%m-%dT%H:%M:%S"
    }
    
    # Logging Handlers
    console_handler = {
        "class"     : "logging.StreamHandler",
        "formatter" : "standard",
        "level"     : console_level
    }
    file_handler = {
        "class"     : "logging.FileHandler",
        "formatter" : "standard",
        "level"     : "INFO",
        "filename"  : cwd + "/motes.log",
		"mode"      : "w",
		"encoding"  : "utf-8"
    }
    
    # MOTES logger
    motes_logger = {
        "handlers"   : ["console", "file"],
        "level"      : "INFO",
         "propagate" : False
    }
	
	# Define Logging Schema
    logging_schema = {
        "disable_existing_loggers" : False,
        "version"    : 1,
        "formatters" : {"standard" : motes_formatter},
        "handlers"   : {"console": console_handler, "file" : file_handler},
        "loggers"    : {"motes" : motes_logger}
    }

    return logging_schema

