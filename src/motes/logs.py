#!/usr/bin/env python3

"""
	logs.py

	Copyright (C) 2025 Tom Seccull & Dominik Kiersz
	
	This module is part of the MOTES package hosted at 
	https://github.com/tseccull/motes
	https://doi.org/####################################################
	
	If used, please cite the MOTES DOI above.
	
	This script is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.

	Last updated - 2025-03-07

	Description---------------------------------------------------------
	logs.py contains a function defining the formatting for MOTES logs.
"""


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
        "format" : "%(asctime)s | %(levelname)8s | %(module)11s | %(funcName)18s | %(message)s",
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

