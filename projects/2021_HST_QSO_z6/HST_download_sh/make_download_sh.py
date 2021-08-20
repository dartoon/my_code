from astroquery.mast import Observations
import numpy as np
import os,sys

search_radius = 120

## import files ##
inp0    = 'ALMA_z6qso_all.list'
data    = np.loadtxt(inp0,comments='#',dtype={'names':('proj_id','name','ra','dec','PI'),'formats':('S32','S32','f16','f16','S32')})
proj_id = data["proj_id"]
name    = data["name"]
ra      = data["ra"]
dec     = data["dec"]

for i in range(len(name)):
	obs_table = Observations.query_criteria(coordinates="%3.5f %3.5f" %(ra[i],dec[i]), radius="%3.5f arcsec" % search_radius, intentType=["science","SCIENCE"], obs_collection=["HST"])  
	__,uidx = np.unique(obs_table['target_name'],return_index=True)
	target_table = obs_table[uidx]['target_name','s_ra','s_dec','filters','t_exptime','proposal_id','dataURL','obsid']
	if len(uidx) > 0:
		os.system('mkdir download_sh/'+name[i])
		for u in range(len(uidx)):
			data_products = Observations.get_product_list(target_table[u][-1])
			Observations.download_products(data_products, calib_level=[2,3], productType="SCIENCE", curl_flag=True,mrp_only=True,download_dir='download_sh/'+name[i])

