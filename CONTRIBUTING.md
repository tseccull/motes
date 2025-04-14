# How to Contribute
First off, thank you for your interest in contributing to MOTES. Its success,
overall usefulness, and level of adoption by the astronomical community is 
enhanced by the input provided by each contributing developer. 

## Questions & Bug Reports
Please bear in mind that this repo is maintained at an absolute minimal level,
free of charge, by one individual within limited amounts of their spare time.
Bug reports and questions about MOTES are welcome, but I cannot guarantee bug 
fixes or timely responses to queries. The more detailed and specific you can be
about any issues you're having with MOTES, the more likely it is that I'll be
able to assist. Please follow the bug reporting guidelines outlined 
[here](https://opensource.guide/how-to-contribute/#how-to-submit-a-contribution) 
to maximise your odds of success. The only accepted changes, additions, or 
extensions of MOTES and its functionality are outlined in the sections below.
Requests for any other changes to MOTES will most likely be ignored.

## Instrument I/O Modules
The only contributions accepted to MOTES are I/O modules for instruments that
are not yet supported. To do this, fork MOTES, work on your new instrument I/O 
module, test it with your own data from that instrument, and request 
integration of the new module via a GitHub pull request. Examples of working 
modules for supported instruments are located in `src/motes/io/`; feel free to 
use these as guides or templates in the development of your own I/O module. To 
ensure that your instrument I/O module can interface with the rest of MOTES, a 
few standard conventions must be followed.

### Names
- For clarity, the name of the module should have the form `instrument_nameio.py`
-     

### Arguments and Returns
 
### Format 

### Testing
