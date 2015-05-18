#!/usr/bin/python

import os
import tempfile

print 'Building - MATLAB - File'

pathTmp = tempfile.gettempdir() + '/%s'
filename1 = 'tmpRTOptFunction'
filename2 = 'tmpRTOptJacobi'
filename3 = 'tmpRTOptHesse'
template = 'template_BasisGenQDyn.m'
file_to = 'BasisGenQDyn.m'

string_to_replace = [
					['q(1)', 'q(1, :)'],
					['q(2)', 'q(2, :)'],
					['q(3)', 'q(3, :)'],
					['q(4)', 'q(4, :)'],
					['v(1)', 'v(1, :)'],
					['v(2)', 'v(2, :)'],
					['v(3)', 'v(3, :)'],
					['omega(1)', 'omega(1, :)'],
					['omega(2)', 'omega(2, :)'],
					['omega(3)', 'omega(3, :)'],
					['u(1)', 'u(1, :)'],
					['u(2)', 'u(2, :)'],
					['u(3)', 'u(3, :)'],
					['u(4)', 'u(4, :)'],
					['*', '.*'],
					['/', './'],
					['^', u'.^'],
					['ones', 'onesV']
					]

path1 = pathTmp % filename1
path2 = pathTmp % filename2
path3 = pathTmp % filename3

if os.path.exists(path1):
	os.remove(path1)
if os.path.exists(path2):
	os.remove(path2)
if os.path.exists(path3):
	os.remove(path3)

os.system('maple < GenerateJacobiHesse.mpl')

template_file = open(template, 'r')
file1 = open(path1, 'r')
file2 = open(path2, 'r')
file3 = open(path3, 'r')
file_write_to = open(file_to, 'w')


try:

	text = template_file.read()
	text1 = file1.read()
	text2 = file2.read()
	text3 = file3.read()

	for str_pair in string_to_replace:
		text1 =  text1.replace(str_pair[0], str_pair[1])

	for str_pair in string_to_replace:
		text2 =  text2.replace(str_pair[0], str_pair[1])

	for str_pair in string_to_replace:
		text3 =  text3.replace(str_pair[0], str_pair[1])
	
	text =  text.format(text1, text2, text3 )
	
	file_write_to.write(text)

	#for str_pair in string_to_replace:
	#	text =  text.replace(str_pair[0], str_pair[1])
	#print text
finally:
	file1.close()
	file2.close()
	file3.close()
	template_file.close()
	file_write_to.close()

 
