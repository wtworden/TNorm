# TNorm

TNorm is a package for computing the Thurston norm unit ball of finite volume orientable hyperbolic 3-manifolds. Currently, tnorm must be installed in Sage. 

There are two ways to get Tnorm up and running. Below is the easy way, scroll down for the hard way (not recommended).

Installing Tnorm:

Install the latest stable version of Sage_macOS, which is a version of Sage packaged by Marc Culler that makes everything work on a Mac without any headaches (thank you, Marc!):

https://github.com/3-manifolds/Sage_macOS/releases/latest  (scroll to bottom and click on SageMath-9.X.dmg)

After opening the installer, drag SageMath-9-X into the Applications folder, then be sure to install Recommended.pkg by double clicking on it and following the prompts.

Next, open a terminal and open sage:

	$: sage

then install tnorm:

	sage: %pip install tnorm

In the course of installing tnorm, Regina and SnapPy will be installed automatically.


To run the tnorm graphical user interface app:

	$ sage -python -m tnorm.app

To get started:

	sage: W=tnorm.load('m130')
	sage: B=W.norm_ball
	sage: B.vertices()
	[Vertex 0: represented by (1/2)* S_1,2 at (-1), mapped from surface with index 10,
	Vertex 1: represented by (1/2)* S_1,2 at (1), mapped from surface with index 0]
	sage: 

Support for hyperbolic 3-manifolds that are not multi-component links in rational homology 3-spheres has been added very recently, and has not been thoroughly tested yet. If you get any results that don't make sense, please email me at william.worden@rice.edu.

TO DO:

* add feature: determine fiberedness of a hyp 3-mfld (this is a long-term goal)
* remove Sage dependence (this has become low-priority, since installing in Sage is now easy).
* better documentation throughout.



The HARD(ish) way (not recommended, but left here to remind us how difficult things were in the old days!):


1. Install latest Sage:
	
	1a. Download from http://mirrors.mit.edu/sage/osx/intel/index.html. Choose the most recent version, ending in x86_64.dmg.
	
	1b. double click on the downloaded file, then after it unpacks, copy the SageMath folder to your Applications folder

2. Unquarantine sage:
	
	2a. open Terminal (in Applications/Utilities)
	
	2b. In terminal, type (without $ sign):
		 
		 $ xattr -rd com.apple.quarantine /Applications/SageMath
		
	(note: for the above, "/Applications/" is the path to SageMath. This may be different if you put it (for example) in the applications folder at the user lever, i.e., /Users/your_username/Applications/SageMath)

	2c. Open System Preferences —> Security and Privacy —> Privacy, then scroll on the left side down to "Full Disk Access" and select it. On the right side (after clicking the lock and entering password, if needed), click “+” button, then navigate to Applications/SageMath/ and select sage.

3. Make a sage alias:
	
	3a. in Terminal, type 
	
		$ open -e .bash_profile
	
	3b. this will open a file in TextEdit. At the top, add the line 
	
		$ alias sage='/Applications/SageMath/sage'
	
	(change the path as needed if your SageMath is in a different location)

4. Fix Sage’s SSL (this is now easy thanks to Nathan Dunfield):

	4a. Download https://github.com/3-manifolds/fix_mac_sage/releases/download/v1.0_as_released/fix_mac_sage9.tar.gz , then double click on it to unpack it.
	
	4b. In Terminal, type: 
	
		$ cd ~/Downloads; sage -python -m fix_mac_sage9.fix

5. Install tnorm (snappy and sageRegina will be installed automatically, if they are not already installed):
	
	5a: In Terminal, type: 
	
		$ sage -python -m pip install tnorm






