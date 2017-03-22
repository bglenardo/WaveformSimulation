from subprocess import call

#with open('tt_means.txt') as tt_m:
#   tt_m_list = [x.strip('\n') for x in tt_m.readlines()]

#with open('tt_sigs.txt') as tt_s:
#   tt_s_list = [x.strip('\n') for x in tt_s.readlines()]

for i in range(0,10):
  call(['./LineupResolutionSingleBin',str(i*10.)])

