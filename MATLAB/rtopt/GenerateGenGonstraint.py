#!/usr/bin/python

import os
import tempfile

print 'Building - MATLAB - File'

pathTmp = tempfile.gettempdir() + '/%s'
filename1 = 'tmpRTConstraintsFunction'
filename2 = 'tmpRTConstraintsJacobi'
filename3 = 'tmpRTConstraintsHesse'
filename4 = 'tmpRTConstraintsJacobiCount'
filename5 = 'tmpRTConstraintsCount'
template = 'template_GenConstraints.m'
file_to = 'GenConstraints.m'

string_to_replace = [
					['q(1)', 'q(1, timestep)'],
					['q(2)', 'q(2, timestep)'],
					['q(3)', 'q(3, timestep)'],
					['q(4)', 'q(4, timestep)'],
					['v(1)', 'v(1, timestep)'],
					['v(2)', 'v(2, timestep)'],
					['v(3)', 'v(3, timestep)'],
					['omega(1)', 'omega(1, timestep)'],
					['omega(2)', 'omega(2, timestep)'],
					['omega(3)', 'omega(3, timestep)'],
					['u(1)', 'u(1, timestep)'],
					['u(2)', 'u(2, timestep)'],
					['u(3)', 'u(3, timestep)'],
					['u(4)', 'u(4, timestep)'],
					['*', '.*'],
					['/', './'],
					['^', u'.^'],
					['ones', 'onesV']
					]

path1 = pathTmp % filename1
path2 = pathTmp % filename2
path3 = pathTmp % filename3
path4 = pathTmp % filename4
path5 = pathTmp % filename5

if os.path.exists(path1):
	os.remove(path1)
if os.path.exists(path2):
	os.remove(path2)
if os.path.exists(path3):
	os.remove(path3)
if os.path.exists(path4):
	os.remove(path4)
if os.path.exists(path5):
	os.remove(path5)

os.system('maple < GenerateConstraints.mpl')

template_file = open(template, 'r')
file1 = open(path1, 'r')
file2 = open(path2, 'r')
file3 = open(path3, 'r')
file4 = open(path4, 'r')
file5 = open(path5, 'r')
file_write_to = open(file_to, 'w')


try:

	text = template_file.read()
	text1 = file1.read()
	text2 = file2.read()
	text3 = file3.read()
	text4 = file4.read()
	text5 = file5.read()

	for str_pair in string_to_replace:
		text1 =  text1.replace(str_pair[0], str_pair[1])

	for str_pair in string_to_replace:
		text2 =  text2.replace(str_pair[0], str_pair[1])

	for str_pair in string_to_replace:
		text3 =  text3.replace(str_pair[0], str_pair[1])

	for str_pair in string_to_replace:
		text4 =  text4.replace(str_pair[0], str_pair[1])

	for str_pair in string_to_replace:
		text5 =  text5.replace(str_pair[0], str_pair[1])
	
	text =  text.format(text1, text2, text3, text4, text5)
	
	file_write_to.write(text)

	#for str_pair in string_to_replace:
	#	text =  text.replace(str_pair[0], str_pair[1])
	#print text
finally:
	file1.close()
	file2.close()
	file3.close()
	file4.close()
	file5.close()
	template_file.close()
	file_write_to.close()

 
