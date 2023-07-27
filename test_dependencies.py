from arnie.mfe import mfe
from arnie.utils import load_package_locations
from arnie.free_energy import free_energy
from DegScore import DegScore
from RiboGraphViz import RGV

seq = 'CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG'

# Arnie mfe utility calls LinearFold for structure. Arnie free_energy utility calls LinearPartition for free energy.

print('Testing Ribotree Dependencies')
print('Last updated June 20, 2021 H Wayment-Steele')
print('\nExample sequence: Hammerhead Ribozyme\n%s\n'%seq)

packages = load_package_locations(DEBUG=True)

print('\n#### Looking for Arnie dependencies ####\n')

necessary_list = ['vienna_2','contrafold_2','eternafold','linearfold','linearpartition']

for pkg in necessary_list:
	if pkg in packages.keys():
		print('Found %s with path %s' % (pkg, packages[pkg]))
	else:
		raise RuntimeError('Could not find package %s' % pkg)

print('\n#### Testing Energy Models ####\n')

count=0

print('Vienna RNAfold:')
print(mfe(seq))
dG_vienna = free_energy(seq)
print(dG_vienna)
if dG_vienna == -15.92:
	print('Vienna dG matches')
else:
	print(dG_vienna, -15.92)
	print('Vienna did not match')
	count += 1


print('\nLinearFold-V:')
print(mfe(seq, linear=True))
dG_linearfoldv = free_energy(seq, linear=True)
if dG_linearfoldv == -15.92:
	print('LinearFold-V dG matches')
else:
	print(dG_linearfoldv, -15.92)
	print('LinearFold-V did not match')
	count +=1

print('\nCONTRAfold (v2_02):')
print(mfe(seq, package='contrafold', viterbi=True)) # setting viterbi=True because CONTRAfold default is MEA structure, not MFE structure
dG_contrafold = free_energy(seq, package='contrafold')
print(dG_contrafold)
if dG_contrafold == -6.87394:
	print('CONTRAfold matches')
else:
	print(dG_contrafold, -6.87394)
	print('CONTRAfold does not match value from CONTRAfold v202 or CONTRAfold-se.')
	count +=1

print('\nLinearFold-C:')
print(mfe(seq, linear=True, package='contrafold'))
dG_linearfoldc= free_energy(seq, linear=True, package='contrafold')
print(dG_linearfoldc)
if dG_linearfoldc==-6.77346:
	print('LinearFold-C matches')
else:
	print(dG_linearfoldc, -6.77346)
	print('LinearFold-C did not match.')
	count +=1

print('\nEternaFold:')
print(mfe(seq, package='eternafold', viterbi=True)) # setting viterbi=True because CONTRAfold default is MEA structure, not MFE structure
dG_eternafold = free_energy(seq, package='eternafold')
print(dG_eternafold)
if dG_eternafold==-13.7489:
	print('EternaFold matches')
else:
	print(dG_eternafold, -13.7489)
	print('EternaFold dG did not match. Is your Arnie eternafold path set to `/path/to/EternaFold/src` ?')
	count +=1

print('\nLinearFold-E:')
print(mfe(seq, linear=True, package='eternafold'))
dG_linearfolde = free_energy(seq, linear=True, package='eternafold')
print(dG_linearfolde)
if dG_linearfolde==-13.28986:
	print('LinearFold-E matches')
else:
	print(dG_eternafold, -13.28986)
	print('LinearFold-E dG did not match.')
	count +=1

print('\n#### Testing DegScore ####\n')

mdl = DegScore(seq, mask_U = False)
print('DegScore for example sequence, mask_U = False: %.3f' % mdl.degscore)
if mdl.degscore == 20.420000000000005:
	print('DegScore matches')
else:
	print(mdl.degscore, 20.420000000000005)
	print('DegScore did not match.')
	count +=1

mdl = DegScore(seq, mask_U=True)
print('DegScore for example sequence, mask_U = True: %.3f' % mdl.degscore)
if mdl.degscore == 16.268000000000004:
	print('DegScore matches')
else:
	print(mdl.degscore, 16.268000000000004)
	print('DegScore did not match.')
	count +=1


print('\n#### Testing RiboGraphViz ####\n')

mdl = RGV(mfe(seq))

print('Building RiboGraphViz object based on mfe structure')
print(mfe(seq))

mdl.run_structure_properties()

if mdl.n_hairpins == 2:
	print('RiboGraphViz hairpin count matches (2)')
else:
	print('RiboGraphViz discrepancy', mdl.n_hairpins, 2)
	count +=1

if count==0:
	print('Congrats! Everything matches.')

else:
	print('%d discrepancies encountered -- see readout above.' % count)
