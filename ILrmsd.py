#该脚本用于计算CABSdock top 10 model的IL_RMSD
#先通过对比complex bindingsites residues和reference_alignment_target.csv，得到complex bindingsites residues与top10 model中氨基酸的对应列表（unbound receptor中缺失的bindingsites residues则不被纳入）
#再计算列表中的top10 model中对应氨基酸在top10 model.pdb成功预测的比例

import pandas as pd
import numpy as np
import os,re,glob
import MDAnalysis as mda
from MDAnalysis.analysis import rms,align


with open('./IL-RMSD.txt',mode='w',encoding='utf-8') as g:
	print('Task\ttop_n\tIL-RMSD',file=g)
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
			
			#提取complex receptor中的interface residues，以及其在top10 model中的对应residues
			#complex中receptor interface residues读取
			interface_resids=[]
			if len(complex_chain_rec)>1:
				for chain in complex_chain_rec:
					df_complex=pd.read_csv('/mnt/hgfs/share/{}/interface_rec_chain{}.csv'.format(task,chain))
					for residue in df_complex.iloc[:,0]:
						if isinstance(residue,str):
							residue=residue.split(' ')
							resid=residue[1]
							interface_resid=chain+':'+resid
							interface_resids.append(interface_resid)
						else: break
			else:
				df_complex=pd.read_csv('/mnt/hgfs/share/{}/interface_rec.csv'.format(task))
				for residue in df_complex.iloc[:,0]:
					if isinstance(residue,str):
						residue=residue.split(' ')
						resid=residue[1]
						interface_resid=complex_chain_rec+':'+resid
						interface_resids.append(interface_resid)
					else: break
			print(interface_resids)
			print(len(interface_resids))

			#CABSdock alignment读取
			df_align=pd.read_csv('/mnt/hgfs/share/{}/reference_alignment_target.csv'.format(task))
			df_align=np.array(df_align)
			reference_residues=[]
			template_residues=[]
			align_residues=[]
			for i in range(df_align.size):
				aligns=str(df_align[i,0])
				aligns=aligns.split('\t')
				reference_residue=aligns[0]
				template_residue=aligns[1]
				reference_residue=reference_residue[:-2]
				template_residue=template_residue[:-2]
				if reference_residue in interface_resids:
					reference_residues.append(reference_residue)
					template_residues.append(template_residue)
					align_residues.append([reference_residue,template_residue])
			print(len(reference_residues))
			print(len(template_residues))
			print(reference_residues)
			print(template_residues)

			#将提取得到的residues写为MDAnalysis中的select语句
			interface_com=''
			if len(complex_chain_rec)>1:
				for chain in complex_chain_rec:
					interface_com_singlechain='segid {} and resid'.format(chain)
					for interface_com_residue in reference_residues:
						if (interface_com_residue[0] == chain):
							interface_com_resid=interface_com_residue[2:]
							interface_com_singlechain=interface_com_singlechain+' '+str(interface_com_resid)
						else:continue
					interface_com_singlechain='('+interface_com_singlechain+')'
					interface_com=interface_com+interface_com_singlechain+' or '
				interface_com=interface_com[:-4]
			else:
				interface_com_singlechain='segid {} and resid'.format(complex_chain_rec)
				for interface_com_residue in reference_residues:
					if (interface_com_residue[0] == complex_chain_rec):
						interface_com_resid=interface_com_residue[2:]
						interface_com_singlechain=interface_com_singlechain+' '+str(interface_com_resid)
					else:continue
				interface_com_singlechain='('+interface_com_singlechain+')'
				interface_com=interface_com_singlechain
			print(interface_com)
			
			interface_model=''
			if len(model_chain_rec)>1:
				for chain in model_chain_rec:
					interface_mod_singlechain='segid {} and resid'.format(chain)
					for interface_mod_residue in template_residues:
						if (interface_mod_residue[0] == chain):
							interface_mod_resid=interface_mod_residue[2:]
							interface_mod_singlechain=interface_mod_singlechain+' '+str(interface_mod_resid)
						else:continue
					interface_mod_singlechain='('+interface_mod_singlechain+')'
					interface_model=interface_model+interface_mod_singlechain+' or '
				interface_model=interface_model[:-4]
			else:
				interface_mod_singlechain='segid {} and resid'.format(model_chain_rec)
				for interface_mod_residue in template_residues:
					if (interface_mod_residue[0] == model_chain_rec):
						interface_mod_resid=interface_mod_residue[2:]
						interface_mod_singlechain=interface_mod_singlechain+' '+str(interface_mod_resid)
					else:continue
				interface_mod_singlechain='('+interface_mod_singlechain+')'
				interface_model=interface_mod_singlechain
			print(interface_model)

			#提取complex peptide中的interface residues，以及其在top10 model中的对应residues
			#complex中peptide interface residues读取
			interface_pep_resids=[]
			df_complex=pd.read_csv('/mnt/hgfs/share/{}/interface_pep.csv'.format(task))
			for residue in df_complex.iloc[:,0]:
				if isinstance(residue,str):
					residue=residue.split(' ')
					resid=residue[1]
					interface_pep_resid=complex_chain_pep+':'+resid
					interface_pep_resids.append(interface_pep_resid)
				else: break
			print(interface_pep_resids)
			print(len(interface_pep_resids))

			#CABSdock alignment读取
			df_align=pd.read_csv('/mnt/hgfs/share/{}/reference_alignment_{}.csv'.format(task,model_chain_pep))
			df_align=np.array(df_align)
			reference_pep_residues=[]
			template_pep_residues=[]
			align_pep_residues=[]
			for i in range(df_align.size):
				aligns=str(df_align[i,0])
				aligns=aligns.split('\t')
				reference_pep_residue=aligns[0]
				template_pep_residue=aligns[1]
				reference_pep_residue=reference_pep_residue[:-2]
				template_pep_residue=template_pep_residue[:-2]
				if reference_pep_residue in interface_pep_resids:
					reference_pep_residues.append(reference_pep_residue)
					template_pep_residues.append(template_pep_residue)
					align_pep_residues.append([reference_pep_residue,template_pep_residue])
			print(len(reference_pep_residues))
			print(len(template_pep_residues))
			print(reference_pep_residues)
			print(template_pep_residues)
			
			#将提取得到的residues写为MDAnalysis中的select语句
			interface_pep_com=''
			interface_pep_com_singlechain='segid {} and resid'.format(complex_chain_pep)
			for interface_pep_com_residue in reference_pep_residues:
				if (interface_pep_com_residue[0] == complex_chain_pep):
					interface_pep_com_resid=interface_pep_com_residue[2:]
					interface_pep_com_singlechain=interface_pep_com_singlechain+' '+str(interface_pep_com_resid)
				else:continue
			interface_pep_com_singlechain='('+interface_pep_com_singlechain+')'
			interface_pep_com=interface_pep_com_singlechain
			print(interface_pep_com)
			
			interface_pep_model=''
			interface_pep_mod_singlechain='segid {} and resid'.format(model_chain_pep)
			for interface_pep_mod_residue in template_pep_residues:
				if (interface_pep_mod_residue[0] == model_chain_pep):
					interface_pep_mod_resid=interface_pep_mod_residue[2:]
					interface_pep_mod_singlechain=interface_pep_mod_singlechain+' '+str(interface_pep_mod_resid)
				else:continue
			interface_pep_mod_singlechain='('+interface_pep_mod_singlechain+')'
			interface_pep_model=interface_pep_mod_singlechain
			print(interface_pep_model)
			
			#sel语句生成后，开始进行rmsd计算
			for i in range(10):
				com=mda.Universe('/mnt/hgfs/share/{}/complex.pdb'.format(task))
				com_backbone=com.select_atoms('backbone')
				model=mda.Universe('/mnt/hgfs/share/{}/model_{}.pdb'.format(task,i))
				model_backbone=model.select_atoms('backbone')
				alignment=align.alignto(com_backbone,model_backbone,select={'mobile':interface_com,'reference':interface_model})
				print(alignment)
				com_pep=com_backbone.select_atoms(interface_pep_com)
				mod_pep=model_backbone.select_atoms(interface_pep_model)
				rmsd=rms.rmsd(com_pep.positions,mod_pep.positions)
				print(task+'\t'+str(i)+'\t'+str(rmsd),file=g)