README.txt


INTRODUCTION:
    Genome Relative Abundance using Mixture Model thoery (GRAMMy).
    Currently the package works for Linux (Ubuntu) platforms.
    It might work for Windows with Cygwin and Mac with Xcode and Mac ports (not tested).
    It requires the Standard C++ Library with Technical Draft One.
    More current information of this package is available @
	http://meta.usc.edu/softs/grammy

DEPENDENCIES:
    C++ build environment
	e.g. build-essential and libstdc++6 in Ubuntu and Xcode in Mac 
    Python(>=2.7) with developemnt and setuptools
        download @ http://www.python.org/
    Numpy(>=1.0)
        download @ http://www.scipy.org/
    Scipy(>=0.6)
        download @ http://www.scipy.org/
    Biopython(>=1.0)
        download @ http://biopython.org/
    
    For setting up the dependencies, users may refer to the author's development document @
	http://dl.dropbox.com/u/35182955/Ubuntu_development_environment.html

FILES:
    LICENSE.txt
    README.txt
    MANIFEST.in
    setup.py
    grammy/*
    test/*

INSTALL:

    [VirtualBox]
    The procedure is similar to QIIME VirtualBox install, 
	see http://qiime.org/install/virtual_box.html.

    1. Download and install the VirtualBox (VB) version for your machine,
	at http://www.virtualbox.org

    2. Download the SunLab Virtual Box, 
	at http://meta.usc.edu/softs/vbox/SunLab.vdi.tgz
	This file is large so it may take
	between a few minutes and a few hours depending on your Internet connection
	speed. You will need to unzip this file, which you can typically do by
	double-clicking on it.

    3. Create a new virtual machine:
        Launch VirtualBox, and create a new machine (press the New button).
        A new window will show up. Click ‘Next’.

        In this screen type SunLab as the name for the virtual machine. Then
	select Linux as the Operating System, and Ubuntu as the version.
	Click Next.

        Select the amount of RAM (memory). You will need at least 512MB, but
	the best option is based on your machine. After selecting the amount of RAM,
	click Next.

        Select “Use existing hard drive”, and click the folder icon next to
	the selector (it has a green up arrow). In the new window click ‘Add’, and
	locate the virtual hard drive that was downloaded in step 2. Click Select and
	then click Next.

        In the new window click Finish.
	
    4. Double click on the new virtual machine created – it will be called SunLab
	– to boot it for the first time. The default username and password is:
	user.

    5. Review any messages that are shown, and select whatever options are best
	for you.

    [Prerequisites]

    Please fullfill the prerequisites of C++, Python (with development and setuptools),
    numpy, scipy and biopython as described in README.txt before installing GRAMMy.

    [Linux] (Ubuntu)

    Download the latest release of GRAMMy from http://meta.usc.edu/softs/grammy.
    Follow standard python module setup to install:
        $tar -zxvf charade-grammy-release.tar.gz
        $cd charade-grammy-release
        $python setup.py install
        $cd test
        $. test.sh

    [Development]

    GRAMMy is open source and the version controled repository is @:
    	https://bitbucket.org/charade/grammy.
    Use mercurial tools (http://mercurial.selenic.com) to download a local copy:
        $hg clone ssh://hg@bitbucket.org/charade/grammy grammy-tip

    Follow standard python module setup to install:
        $cd grammy-tip
        $python setup.py install

EXECUTABLES:
    grammy_gdt
    grammy_rdt
    grammy_pre
    grammy_em
    grammy_post

USAGE:
    (i) Above executables will be available from your python scripts directory.
    	Use '-h' to read individual script usage.
    (ii) A pre-structured genome and taxon directory is required for grammy_gdt.py;
    	An example is available as the grefs.tgz file.
    	The file can be downloaded from http://meta.usc.edu/softs/grammy/grefs.tgz
    (iii) A simple test example is available at 'test/test.sh' and explained with in.

