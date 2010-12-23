# toolbox setup
def settings():
	import platform
	
	sw_base = "/data/BI/Toolbox/software/"
	
	if platform.system() == "Darwin":
		if platform.architecture()[0] == "64bit":
			#MAC 64 bit
			fsl_path  = sw_base + "fsl_maci_4.1.5/bin/"
			ants_path = sw_base + "ants_maci_svn460/bin/"
			afni_path = sw_base + "afni_maci_os10.5_64bit_v1710/"
			c3d_path  = sw_base + "c3d_mac_maci_0.8.0/bin/"
		elif platform.architecture()[0] == "32bit":
			#MAC 32 bit
			fsl_path  = sw_base + "fsl_maci_4.1.5/bin/"
			ants_path = sw_base + "ants_maci_svn460/bin/"
			afni_path = sw_base + "afni_maci_os10.5_32bit_v1710/"
			c3d_path  = sw_base + "c3d_mac_maci_0.8.0/bin/"
		else:
			print("Error: running on unsupported architecture " + platform.system())
			quit()
	
	elif platform.system() == "Linux":
		if platform.architecture()[0] == "64bit":
			#LINUX 64 bit
			fsl_path  = sw_base + "fsl_lnx64_v4.1.4/bin/"
			ants_path = sw_base + "ants_lnx64_svn460/bin/"
			afni_path = sw_base + "afni_lnx_2009_12_31_1431/"
			c3d_path  = sw_base + "c3d_lnx64_0.8.2/bin/"
		elif platform.architecture()[0] == "32bit":
			#LINUX 32 bit
			fsl_path  = sw_base + "fsl_lnx_4.1.4/bin/"
			ants_path = sw_base + ""
			afni_path = sw_base + "afni_lnx_2009_12_31_1431/"
			c3d_path  = sw_base + ""
		else:
			print("Error: running on unsupported architecture " + platform.system())
			quit()		
	return [fsl_path,ants_path,afni_path,c3d_path]