pro benchmark06, windx=WINDX, oname=ONAME
;
;plots contours of searchlight along a direction
;
;------------------read all information from hdf5-file------------------
;
fname='benchmark06.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'nd_fine')
         nd_fine=h5a_read(att_id)
         nd_fine=nd_fine(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nd_method01')
         nd_method01=h5a_read(att_id)
         nd_method01=nd_method01(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nd_method02')
         nd_method02=h5a_read(att_id)
         nd_method02=nd_method02(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nd_method03')
         nd_method03=h5a_read(att_id)
         nd_method03=nd_method03(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nd_method04')
         nd_method04=h5a_read(att_id)
         nd_method04=nd_method04(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndc_fine')
         ndc_fine=h5a_read(att_id)
         ndc_fine=ndc_fine(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndc_method01')
         ndc_method01=h5a_read(att_id)
         ndc_method01=ndc_method01(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndc_method02')
         ndc_method02=h5a_read(att_id)
         ndc_method02=ndc_method02(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndc_method03')
         ndc_method03=h5a_read(att_id)
         ndc_method03=ndc_method03(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndc_method04')
         ndc_method04=h5a_read(att_id)
         ndc_method04=ndc_method04(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method00')
      dset_id=h5d_open(group_id, 's1d_fine')
         s1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1d_fine')
         vel1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_fine')
         opalbar1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1d_fine')
         vth1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1d_fine')
         sline1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1d_fine')
         int1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1d_fine')
         opal1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1d_fine')
         tau1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1d_fine')
         profile1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1d_fine')
         xcmf1d_fine=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method01')
      dset_id=h5d_open(group_id, 's1d_method01')
         s1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1d_method01')
         vel1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_method01')
         opalbar1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1d_method01')
         vth1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1d_method01')
         sline1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1d_method01')
         int1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1d_method01')
         opal1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1d_method01')
         tau1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1d_method01')
         profile1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1d_method01')
         xcmf1d_method01=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method02')
      dset_id=h5d_open(group_id, 's1d_method02')
         s1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1d_method02')
         vel1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_method02')
         opalbar1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1d_method02')
         vth1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1d_method02')
         sline1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1d_method02')
         int1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1d_method02')
         opal1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1d_method02')
         tau1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1d_method02')
         profile1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1d_method02')
         xcmf1d_method02=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method03')
      dset_id=h5d_open(group_id, 's1d_method03')
         s1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1d_method03')
         vel1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_method03')
         opalbar1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1d_method03')
         vth1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1d_method03')
         sline1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1d_method03')
         int1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1d_method03')
         opal1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1d_method03')
         tau1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1d_method03')
         profile1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1d_method03')
         xcmf1d_method03=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method04')
      dset_id=h5d_open(group_id, 's1d_method04')
         s1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1d_method04')
         vel1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_method04')
         opalbar1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1d_method04')
         vth1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1d_method04')
         sline1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1d_method04')
         int1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1d_method04')
         opal1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1d_method04')
         tau1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1d_method04')
         profile1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1d_method04')
         xcmf1d_method04=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method00c')
      dset_id=h5d_open(group_id, 's1dc_fine')
         s1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1dc_fine')
         vel1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1dc_fine')
         opalbar1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1dc_fine')
         vth1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1dc_fine')
         sline1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont1dc_fine')
         scont1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'stot1dc_fine')
         stot1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1dc_fine')
         int1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1dc_fine')
         opal1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1dc_fine')
         opac1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opatot1dc_fine')
         opatot1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1dc_fine')
         tau1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1dc_fine')
         profile1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1dc_fine')
         xcmf1dc_fine=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method01c')
      dset_id=h5d_open(group_id, 's1dc_method01')
         s1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1dc_method01')
         vel1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1dc_method01')
         opalbar1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1dc_method01')
         vth1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1dc_method01')
         sline1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont1dc_method01')
         scont1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'stot1dc_method01')
         stot1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1dc_method01')
         int1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1dc_method01')
         opal1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1dc_method01')
         opac1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opatot1dc_method01')
         opatot1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1dc_method01')
         tau1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1dc_method01')
         profile1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1dc_method01')
         xcmf1dc_method01=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method02c')
      dset_id=h5d_open(group_id, 's1dc_method02')
         s1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1dc_method02')
         vel1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1dc_method02')
         opalbar1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1dc_method02')
         vth1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1dc_method02')
         sline1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont1dc_method02')
         scont1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'stot1dc_method02')
         stot1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1dc_method02')
         int1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1dc_method02')
         opal1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1dc_method02')
         opac1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opatot1dc_method02')
         opatot1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1dc_method02')
         tau1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1dc_method02')
         profile1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1dc_method02')
         xcmf1dc_method02=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method03c')
      dset_id=h5d_open(group_id, 's1dc_method03')
         s1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1dc_method03')
         vel1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1dc_method03')
         opalbar1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1dc_method03')
         vth1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1dc_method03')
         sline1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont1dc_method03')
         scont1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'stot1dc_method03')
         stot1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1dc_method03')
         int1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1dc_method03')
         opal1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1dc_method03')
         opac1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opatot1dc_method03')
         opatot1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1dc_method03')
         tau1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1dc_method03')
         profile1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1dc_method03')
         xcmf1dc_method03=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'method04c')
      dset_id=h5d_open(group_id, 's1dc_method04')
         s1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vel1dc_method04')
         vel1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1dc_method04')
         opalbar1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth1dc_method04')
         vth1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1dc_method04')
         sline1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont1dc_method04')
         scont1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'stot1dc_method04')
         stot1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int1dc_method04')
         int1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opal1dc_method04')
         opal1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1dc_method04')
         opac1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opatot1dc_method04')
         opatot1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau1dc_method04')
         tau1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'profile1dc_method04')
         profile1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xcmf1dc_method04')
         xcmf1dc_method04=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;----------------------------title-strings------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
xtitlestr=textoidl('s')
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark06.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $; xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx, xsize=1600, ysize=1600
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
;-----------------------------------------------------------------------
;
xlim=[min(s1d_fine),max(s1d_fine)]
xtitlestr=textoidl('s')
;
!p.multi=[20,5,4]
plot, s1dc_fine, opac1dc_fine, xrange=xlim, yrange=[min(opac1dc_fine),max(opac1dc_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('\chi_c')
oplot, s1dc_method01, opac1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, opac1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, opac1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, opac1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[19,5,4]
plot, s1d_fine, opalbar1d_fine, xrange=xlim, yrange=[min(opalbar1d_fine),max(opalbar1d_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('<\chi_L>')
oplot, s1d_method01, opalbar1d_method01, line=2, color=ci_blue
oplot, s1d_method02, opalbar1d_method02, line=2, color=ci_red
oplot, s1d_method03, opalbar1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, opalbar1d_method04, line=2, color=ci_green
oplot, s1dc_fine, opalbar1dc_fine, thick=2.
oplot, s1dc_method01, opalbar1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, opalbar1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, opalbar1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, opalbar1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[18,5,4]
ylim=[min([opal1d_fine,opal1d_method02,opal1d_method04]),max([opal1d_fine,opal1d_method02,opal1d_method04])]
plot, s1d_fine, opal1d_fine, xrange=xlim, yrange=ylim, charsize=2., xtitle=xtitlestr, ytitle=textoidl('\chi_L')
oplot, s1d_method01, opal1d_method01, line=2, color=ci_blue
oplot, s1d_method02, opal1d_method02, line=2, color=ci_red
oplot, s1d_method03, opal1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, opal1d_method04, line=2, color=ci_green
oplot, s1dc_fine, opal1dc_fine, thick=2.
oplot, s1dc_method01, opal1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, opal1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, opal1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, opal1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[17,5,4]
ylim=[min([opatot1dc_fine,opatot1dc_method02]),max([opatot1dc_fine,opatot1dc_method02])]
plot, s1dc_fine, opatot1dc_fine, xrange=xlim, yrange=ylim, charsize=2., xtitle=xtitlestr, ytitle=textoidl('\chi_{tot}'), thick=2.
oplot, s1dc_method01, opatot1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, opatot1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, opatot1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, opatot1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[16,5,4]
lstr0='without cont,     fine grid, linear'
lstr1='without cont, no refinement, linear'
lstr2='without cont,    refinement, linear'
lstr3='without cont, no refinement, bezier'
lstr4='without cont,    refinement, bezier'
lstr0c='with cont,     fine grid, linear'
lstr1c='with cont, no refinement, linear'
lstr2c='with cont,    refinement, linear'
lstr3c='with cont, no refinement, bezier'
lstr4c='with cont,    refinement, bezier'
plot, [0,0], [0,0], charsize=2., xticks=1, yticks=1, xtickformat='(a1)', ytickformat='(a1)'
if(keyword_set(oname)) then begin
   legend, [lstr0, lstr1, lstr2, lstr3, lstr4, lstr0c, lstr1c, lstr2c, lstr3c, lstr4c],  $
           psym=[0,0,0,0,0,0,0,0,0,0], $
           linestyle=[0,2,2,2,2,0,2,2,2,2], $
           thick=[0,0,0,0,0,2,2,2,2,2], $
           color=[ci_black,ci_blue, ci_red, ci_cyan, ci_green, ci_black,ci_blue, ci_red, ci_cyan, ci_green], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr0, lstr1, lstr2, lstr3, lstr4, lstr0c, lstr1c, lstr2c, lstr3c, lstr4c],  $
           psym=[0,0,0,0,0,0,0,0,0,0], $
           linestyle=[0,2,2,2,2,0,2,2,2,2], $
           thick=[0,0,0,0,0,2,2,2,2,2], $
           color=[ci_white,ci_blue, ci_red, ci_cyan, ci_green, ci_white, ci_blue, ci_red, ci_cyan, ci_green], $
           textcolor=ci_white, $
           /right_legend
endelse
;
!p.multi=[15,5,4]
plot, s1dc_fine, scont1dc_fine, xrange=xlim, yrange=[min(scont1dc_fine),max(scont1dc_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('S_c'), thick=2.
oplot, s1dc_method01, scont1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, scont1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, scont1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, scont1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[14,5,4]
plot, s1d_fine, sline1d_fine, xrange=xlim, yrange=[min(sline1d_fine),max(sline1d_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('S_L'), /xs, /ys
oplot, s1d_method01, sline1d_method01, line=2, color=ci_blue
oplot, s1d_method02, sline1d_method02, line=2, color=ci_red
oplot, s1d_method03, sline1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, sline1d_method04, line=2, color=ci_green
oplot, s1dc_fine, sline1dc_fine, thick=2.
oplot, s1dc_method01, sline1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, sline1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, sline1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, sline1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[13,5,4]
plot, s1dc_fine, stot1dc_fine, xrange=xlim, yrange=[min(stot1dc_fine),max(stot1dc_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('S_{tot}'), thick=2.
oplot, s1dc_method01, stot1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, stot1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, stot1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, stot1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[12,5,4]
!p.multi=[11,5,4]
;
!p.multi=[10,5,4]
plot, s1d_fine, vth1d_fine, xrange=xlim, yrange=[min(vth1d_fine),max(vth1d_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('v_{th}')
oplot, s1d_method01, vth1d_method01, line=2, color=ci_blue
oplot, s1d_method02, vth1d_method02, line=2, color=ci_red
oplot, s1d_method03, vth1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, vth1d_method04, line=2, color=ci_green
oplot, s1dc_fine, vth1dc_fine, thick=2.
oplot, s1dc_method01, vth1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, vth1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, vth1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, vth1dc_method04, line=2, color=ci_green, thick=2.

!p.multi=[9,5,4]
plot, s1d_fine, vel1d_fine, xrange=xlim, yrange=[min(vel1d_fine),max(vel1d_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('v_s/v_{th}^*')
oplot, s1d_method01, vel1d_method01, line=2, color=ci_blue
oplot, s1d_method02, vel1d_method02, line=2, color=ci_red
oplot, s1d_method03, vel1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, vel1d_method04, line=2, color=ci_green
oplot, s1dc_fine, vel1dc_fine, thick=2.
oplot, s1dc_method01, vel1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, vel1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, vel1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, vel1dc_method04, line=2, color=ci_green, thick=2.

!p.multi=[8,5,4]
plot, s1d_fine, xcmf1d_fine, xrange=xlim, yrange=[min(xcmf1d_fine),max(xcmf1d_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('x_{CMF}')
oplot, s1d_method01, xcmf1d_method01, line=2, color=ci_blue
oplot, s1d_method02, xcmf1d_method02, line=2, color=ci_red
oplot, s1d_method03, xcmf1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, xcmf1d_method04, line=2, color=ci_green
oplot, s1dc_fine, xcmf1dc_fine, thick=2.
oplot, s1dc_method01, xcmf1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, xcmf1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, xcmf1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, xcmf1dc_method04, line=2, color=ci_green, thick=2.

!p.multi=[7,5,4]
plot, s1d_fine, profile1d_fine, xrange=xlim, yrange=[min(profile1d_fine),max(profile1d_fine)], charsize=2., xtitle=xtitlestr, ytitle=textoidl('\Phi')
oplot, s1d_method01, profile1d_method01, line=2, color=ci_blue
oplot, s1d_method02, profile1d_method02, line=2, color=ci_red
oplot, s1d_method03, profile1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, profile1d_method04, line=2, color=ci_green
oplot, s1dc_fine, profile1dc_fine, thick=2.
oplot, s1dc_method01, profile1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, profile1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, profile1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, profile1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[6,5,4]
ylim=[min([int1d_fine,int1d_method01,int1d_method02,int1d_method03,int1d_method04]), $
      max([int1d_fine,int1d_method01,int1d_method02,int1d_method03,int1d_method04])]
plot, s1d_fine, int1d_fine, xrange=xlim, yrange=ylim, charsize=2., xtitle=xtitlestr, ytitle=textoidl('I')
oplot, s1d_method01, int1d_method01, line=2, color=ci_blue
oplot, s1d_method02, int1d_method02, line=2, color=ci_red
oplot, s1d_method03, int1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, int1d_method04, line=2, color=ci_green
oplot, s1dc_fine, int1dc_fine, thick=2.
oplot, s1dc_method01, int1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, int1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, int1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, int1dc_method04, line=2, color=ci_green, thick=2.
;for i=0, nd_method04-1 do begin
;  oplot, [s1d_method04(i),s1d_method04(i)], [min(int1d_fine),max(int1d_fine)], color=ci_green
;endfor

;
!p.multi=[5,5,4]
ylim=[min([tau1d_fine,tau1d_method01,tau1d_method02,tau1d_method03,tau1d_method04]), $
      max([tau1d_fine,tau1d_method01,tau1d_method02,tau1d_method03,tau1d_method04])]
plot, s1d_fine, tau1d_fine, xrange=xlim, yrange=ylim, charsize=2., xtitle=xtitlestr, ytitle=textoidl('\tau')
oplot, s1d_method01, tau1d_method01, line=2, color=ci_blue
oplot, s1d_method02, tau1d_method02, line=2, color=ci_red
oplot, s1d_method03, tau1d_method03, line=2, color=ci_cyan
oplot, s1d_method04, tau1d_method04, line=2, color=ci_green
oplot, s1dc_fine, tau1dc_fine, thick=2.
oplot, s1dc_method01, tau1dc_method01, line=2, color=ci_blue, thick=2.
oplot, s1dc_method02, tau1dc_method02, line=2, color=ci_red, thick=2.
oplot, s1dc_method03, tau1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, s1dc_method04, tau1dc_method04, line=2, color=ci_green, thick=2.
;
!p.multi=[4,5,4]
stot1d_fine=sline1d_fine
stot1d_method01=sline1d_method01
stot1d_method02=sline1d_method02
stot1d_method03=sline1d_method03
stot1d_method04=sline1d_method04
ylim=[min(stot1d_fine),max(stot1d_fine)]
plot, tau1d_fine, sline1d_fine, yrange=ylim, charsize=2., xtitle=textoidl('\tau'), ytitle=textoidl('S_{tot}')
oplot, tau1d_method01, stot1d_method01, line=2, color=ci_blue
oplot, tau1d_method02, stot1d_method02, line=2, color=ci_red
oplot, tau1d_method03, stot1d_method03, line=2, color=ci_cyan
oplot, tau1d_method04, stot1d_method04, line=2, color=ci_green
oplot, tau1dc_fine, stot1dc_fine, thick=2.
oplot, tau1dc_method01, stot1dc_method01, line=2, color=ci_blue, thick=2.
oplot, tau1dc_method02, stot1dc_method02, line=2, color=ci_red, thick=2.
oplot, tau1dc_method03, stot1dc_method03, line=2, color=ci_cyan, thick=2.
oplot, tau1dc_method04, stot1dc_method04, line=2, color=ci_green, thick=2.
;for i=0, nd_method04-1 do begin
;  oplot, [tau1d_method04(i),tau1d_method04(i)], [min(stot1d_fine),max(stot1d_fine)], color=ci_green
;endfor

;
!p.multi=[3,5,4]
plot, xcmf1d_fine, profile1d_fine, yrange=[min(profile1d_fine),max(profile1d_fine)], charsize=2., xtitle=textoidl('x_{CMF}'), ytitle=textoidl('\Phi')
oplot, xcmf1d_method01, profile1d_method01, line=2, color=ci_blue
for i=0, nd_method01-1 do begin
   oplot, [xcmf1d_method01(i),xcmf1d_method01(i)], [min(profile1d_fine),max(profile1d_fine)], color=ci_blue
endfor
oplot, xcmf1d_method02, profile1d_method02, line=2, color=ci_red
for i=0, nd_method02-1 do begin
   oplot, [xcmf1d_method02(i),xcmf1d_method02(i)], [min(profile1d_fine),max(profile1d_fine)], color=ci_red
endfor
oplot, xcmf1d_method03, profile1d_method03, line=2, color=ci_cyan
for i=0, nd_method03-1 do begin
   oplot, [xcmf1d_method03(i),xcmf1d_method03(i)], [min(profile1d_fine),max(profile1d_fine)], color=ci_cyan
endfor
oplot, xcmf1d_method04, profile1d_method04, line=2, color=ci_green
for i=0, nd_method04-1 do begin
   oplot, [xcmf1d_method04(i),xcmf1d_method04(i)], [min(profile1d_fine),max(profile1d_fine)], color=ci_green
endfor
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
!p.multi=0
;
;
end
