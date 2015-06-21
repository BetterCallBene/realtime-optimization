#!/usr/bin/python

import sys
import GenerateScript as GS

MAPLE_PATH = '/Library/Frameworks/Maple.framework/Versions/Current/bin/maple'

def main(template):
	
	files = [
			('tmpRTEqConstraintsFunction', 1),  # 0
			('tmpRTEqConstraintsJacobi',   1),  # 1
			('tmpRTEqConstraintsHesse',    1),  # 2
			('tmpRTEqCountConstraintsJacobi', 1),# 3
			('tmpRTEqCountConstraints', 1),     # 4
			('tmpRTEqCountConstraintsHesse', 1), # 5
			('tmpRTInEqConstraintsFunction', 1), # 6
			('tmpRTInEqConstraintsJacobi', 1),  # 7
			('tmpRTInEqConstraintsHesse', 1),   # 8
			('tmpRTInEqCountConstraintsJacobi', 1), # 9
			('tmpRTInEqCountConstraints', 1),    # 10
			('tmpRTInEqCountConstraintsHesse', 1) # 11
			]
	if template is None:
		template = 'template_GenConstraints.m'

	file_to = 'GenConstraints.m'
	file_mpl = 'GenerateConstraints.mpl'

	strings_to_replace = [
						['r(1)', 'r(1, timestep)'],
						['r(2)', 'r(2, timestep)'],
						['r(3)', 'r(3, timestep)'],
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
	
	GS.generate(file_to, files, file_mpl, template, strings_to_replace)


if __name__ == "__main__":
	template = ''
	if sys.argv[1] is not None:
		template = sys.argv[1]

	main(template)

 
