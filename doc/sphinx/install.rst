.. include:: defs.rst

==============
 Installation
==============

|tenscalc| is available on |github| and requires the installation of 3
toolboxes, which must be installed in the following order:

1. |funpartools|: `<https://github.com/hespanha/funpartools>`_
2. |cmextools|:  `<https://github.com/hespanha/cmextools>`_
3. |tenscalc|:  `<https://github.com/hespanha/tenscalc>`_

Installations instructions for all these tools is available at
|github|.

Issues
======

While most |matlab| scripts are agnostic to the underlying operating
systems (OSs), the use of mex functions depends heavily on the
operating system.

Our goal is to build a toolbox that works across multiple OSs; at
least under OSX, linux, and Microsoft Windows. However, most of our
testing was done under OSX so one should expect some bugs under the
other OSs. Currently, it is fair to say that |tenscalc| has been tested

* fairly extensively under OSX
* lightly under linux
* very lightly under Microsoft Windows

Any help in fixing bugs is greatly appreciated.

  
