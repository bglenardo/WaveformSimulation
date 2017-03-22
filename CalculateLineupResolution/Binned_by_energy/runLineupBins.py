from subprocess import call

with open('dd_means.txt') as dd_m:
   dd_m_list = [x.strip('\n') for x in dd_m.readlines()]

with open('dd_sigs.txt') as dd_s:
   dd_s_list = [x.strip('\n') for x in dd_s.readlines()]

with open('tt_means.txt') as tt_m:
   tt_m_list = [x.strip('\n') for x in tt_m.readlines()]

with open('tt_sigs.txt') as tt_s:
   tt_s_list = [x.strip('\n') for x in tt_s.readlines()]

for i in range(0,len(dd_m_list)):
  call(['./LineupResolutionSingleBin',dd_m_list[i],dd_s_list[i]])

for i in range(0,len(tt_m_list)):
  call(['./LineupResolutionSingleBin',tt_m_list[i],tt_s_list[i]])
