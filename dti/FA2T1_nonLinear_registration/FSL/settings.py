# toolbox setup
import platform

sw_base = "/data/BI/Toolbox/software/"

if platform.system() == "Darwin":
	# Mac settings
	fsl_path  = sw_base + "fsl_maci_4.1.5/bin/"
	ANTS_path = sw_base + "ants_maci_svn460/bin/"
	c3d_path  = sw_base + "c3d_mac_maci_0.8.0/bin/"
	
	if platform.architecture()[0] == "64bit":
		afni_path = sw_base + "afni_maci_os10.5_64bit_v1710/"
	elif platform.architecture()[0] == "32bit":
		afni_path = sw_base + "afni_maci_os10.5_32bit_v1710/"
	else:
		print("Error: running on unsupported architecture " + platform.system())
		quit()
elif platform.system() == "Linux":
	# Linux settings
	fsl_path  = sw_base + "fsl_lnx_4.1.4/bin/"
	ANTS_path = sw_base + "ants_lnx64_svn460/bin/"
	c3d_path  = "No c3d folder for linux"
	afni_path = sw_base + "afni_lnx_2009_12_31_1431/"
else:
	print("Error: running on unsupported architecture " + platform.system())
	quit()