#该脚本用于计算complex.pdb中bindingsites residues被top10 model正确预测的比例
#先通过对比complex bindingsites residues和reference_alignment_target.csv，得到complex bindingsites residues与top10 model中氨基酸的对应列表（unbound receptor中缺失的bindingsites residues则不被纳入）
#再计算列表中的top10 model中对应氨基酸在top10 model.pdb成功预测的比例

import pandas as pd
import numpy as np
import os,re,glob


with open('./bindingsites_predict_proportion.txt',mode='w',encoding='utf-8') as g:
	print('Task\ttop_n\tbindingsite_complex\tbindingsite_alignment\tbindingsite_unboundrec\tcorrect_pre\tproportion_pre',file=g)
	with open('./benchmark_runlist_withouttitle1.txt',encoding='utf-8') as f:
		for line in f:
			code=line.strip('\n')
			item=code.split()
			task=item[0]
			print(task)
			#complex.pdb中peptide和receptor的chainID提取
			complex_chain_rec=item[1]
			complex_chain_pep=item[2]
			model_chain_rec=item[4]
			model_chain_pep=item[3]
			
			#complex中bindingsites residues读取
			bindingsite_resids=[]
			if len(complex_chain_rec)>1:
				for chain in complex_chain_rec:
					df_complex=pd.read_csv('D:/global_docking/benchmark/DP_benchmark/{}/bindingsites_residues_chain{}.csv'.format(task,chain))
					for residue in df_complex.iloc[:,0]:
						if isinstance(residue,str):
							residue=residue.split(' ')
							resid=residue[1]
							bindingsite_resid=chain+':'+resid
							bindingsite_resids.append(bindingsite_resid)
						else: break
			else:
				df_complex=pd.read_csv('D:/global_docking/benchmark/DP_benchmark/{}/bindingsites_residues.csv'.format(task))
				for residue in df_complex.iloc[:,0]:
					if isinstance(residue,str):
						residue=residue.split(' ')
						resid=residue[1]
						bindingsite_resid=complex_chain_rec+':'+resid
						bindingsite_resids.append(bindingsite_resid)
					else: break				   
			print(bindingsite_resids)
			print(len(bindingsite_resids))

			#CABSdock alignment读取
			df_align=pd.read_csv('D:/global_docking/benchmark/CABSdock_PSI_DPbenchmark/'+task+'/'+task+'/output_data/reference_alignment_target.csv')
			df_align=np.array(df_align)
			reference_residues=[]
			template_residues=[]
			align_residues=[]
			for i in range(df_align.size):
				align=str(df_align[i,0])
				align=align.split('\t')
				reference_residue=align[0]
				template_residue=align[1]
				reference_residue=reference_residue[:-2]
				template_residue=template_residue[:-2]
				if reference_residue in bindingsite_resids:
					reference_residues.append(reference_residue)
					template_residues.append(template_residue)
					align_residues.append([reference_residue,template_residue])
			print(len(reference_residues))
			print(len(template_residues))
			print(template_residues)

			if len(model_chain_rec)>1:
				for i in range(10):
					model_bindingsites=[]
					for chain in model_chain_rec:
						df_model=pd.read_csv('D:/global_docking/benchmark/CABSdock_PSI_DPbenchmark/{}/{}/top10_model/model_{}_bindingsites_residues_chain{}.csv'.format(task,task,i,chain))
						for residue in df_model.iloc[:,0]:
							if isinstance(residue,str):
								residue=residue.split(' ')
								resid=residue[1]
								bindingsite_resid=chain+':'+resid
								print(bindingsite_resid)
								if bindingsite_resid in template_residues:
									model_bindingsites.append(bindingsite_resid)
								else:continue
							else:break
					print('*******'+str(i))
					print(model_bindingsites)
					print(len(model_bindingsites))
					print(task+'\t'+str(i)+'\t'+str(len(bindingsite_resids))+'\t'+str(len(reference_residues))+'\t'+str(len(template_residues))+'\t'+str(len(model_bindingsites))+'\t'+str(len(model_bindingsites)/len(reference_residues)),file=g)
			else:
				for i in range(10):
					model_bindingsites=[]
					df_model=pd.read_csv('D:/global_docking/benchmark/CABSdock_PSI_DPbenchmark/{}/{}/top10_model/model_{}_bindingsites_residues.csv'.format(task,task,i))
					for residue in df_model.iloc[:,0]:
						if isinstance(residue,str):
							residue=residue.split(' ')
							resid=residue[1]
							bindingsite_resid=str(model_chain_rec)+':'+resid
							print(bindingsite_resid)
							if bindingsite_resid in template_residues:
								model_bindingsites.append(bindingsite_resid)
							else:continue
						else:break
					print('*******'+str(i))
					print(model_bindingsites)
					print(len(model_bindingsites))
					print(task+'\t'+str(i)+'\t'+str(len(bindingsite_resids))+'\t'+str(len(reference_residues))+'\t'+str(len(template_residues))+'\t'+str(len(model_bindingsites))+'\t'+str(len(model_bindingsites)/len(reference_residues)),file=g)