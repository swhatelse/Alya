require '/tmp/alya.rb'
require 'narray_ffi'
nests = [7]
vector_size=1

k_orig_params = {:vector_length => vector_size, :preprocessor => false, :nests => nests}
k_boast_params = {:vector_length => vector_size, :nests => nests, :unroll => true}

set_fortran_line_length(100)
k_orig = generate_ref_v2(k_orig_params)
file_orig = File::new("orig.f90","w")

if file_orig then
  file_orig.syswrite(k_orig)
else
  puts "Failed to open orig.f90"
end
k_orig.build(:FCFLAGS => "-cpp")

k_boast = generate_boast_implem(k_boast_params)
file_boast = File::new("boast.f90","w")
if file_boast then
  file_boast.syswrite(k_boast)
else
  puts "Failed to open boast.f90"
end
puts k_boast
k_boast.build

kfl_lumped = 1
mnode = 10
ntens = 2
kfl_duatss = 1
fact_duatss = 2
kfl_stabi_nsi = -1 # or 1
fvins_nsi = 2.0
fcons_nsi = 0.5
bemol_nsi = 1.0
kfl_regim_nsi = 3
fvela_nsi = ANArray.float(64,3).random!
kfl_rmom2_nsi = 1
kfl_press_nsi = 1
kfl_p1ve2_nsi = 1
kfl_linea_nsi = 2

kfl_confi_nsi = 1
nbdfp_nsi = 5
kfl_sgsti_nsi = 1
kfl_nota1_nsi = 0
kfl_limit_nsi = 1 # or 2
kfl_penal_nsi = 1
penal_nsi = 1.0
kfl_bubbl_nsi = 1

pnode = 10
pgaus = 10
pevat = 10
ndime = 2
fvins_nsi = 2.0
kfl_sgsti_nsi = 1
nbdfp_nsi = 3

stats_boast = []
stats_ref = []

ANArray.srand(10)
# # 100.times{|i|
i = 1
gpden = ANArray.float(64,vector_size,pgaus).random!
gpvis = ANArray.float(64,vector_size,pgaus).random!
gppor = ANArray.float(64,vector_size,pgaus).random! 
gpsp1 = ANArray.float(64,vector_size,pgaus).random! 
gpsp2 = ANArray.float(64,vector_size,pgaus).random! 
gpvol = ANArray.float(64,vector_size,pgaus).random! 
gpsha = ANArray.float(64,vector_size,pnode,pgaus).random! 
gpcar = ANArray.float(64,vector_size,ndime,mnode,pgaus).random!
gpadv = ANArray.float(64,vector_size,ndime,pgaus).random! 
gpvep = ANArray.float(64,vector_size,ndime,pgaus).random! 
gpprp = ANArray.float(64,vector_size,pgaus).random! 
gpgrp = ANArray.float(64,vector_size,ndime,pgaus).random! 
gprhs = ANArray.float(64,vector_size,ndime,pgaus).random!
gprhc = ANArray.float(64,vector_size,pgaus).random! 
gpvel = ANArray.float(64,vector_size,ndime,pgaus,10).random! #dynamic
gpsgs = ANArray.float(64,vector_size,ndime,pgaus,10).random! #dynamic 
elvel = ANArray.float(64,vector_size,ndime,pnode,10).random! #dynamic
elpre = ANArray.float(64,vector_size,pnode,10).random!
elbub = ANArray.float(64,vector_size).random!

dtinv_loc = ANArray.float(64,vector_size).random!
dtsgs = ANArray.float(64,vector_size).random!
pbubl = ANArray.int(64,vector_size).random!
gpsha_bub = ANArray.float(64,vector_size,pgaus).random!
gpcar_bub = ANArray.float(64,vector_size,ndime,pgaus).random!

wgrgr_ref = ANArray.float(64,vector_size,pnode,pnode,pgaus)
agrau_ref = ANArray.float(64,vector_size,pnode,pgaus)
elauu_ref = ANArray.float(64,vector_size,pnode*ndime,pnode*ndime)
elaup_ref = ANArray.float(64,vector_size,pnode*ndime,pnode)
elapp_ref = ANArray.float(64,vector_size,pnode,pnode)
elapu_ref = ANArray.float(64,vector_size,pnode,pnode*ndime)
elrbu_ref = ANArray.float(64,vector_size,ndime,pnode)
elrbp_ref = ANArray.float(64,vector_size,pnode)

elauq_ref = ANArray.float(vector_size,pnode*ndime,1)
elapq_ref = ANArray.float(vector_size,pnode,1)
elaqu_ref = ANArray.float(vector_size,1,pnode*ndime)
elaqp_ref = ANArray.float(vector_size,1,pnode)
elaqq_ref = ANArray.float(vector_size,1,1)
elrbq_ref = ANArray.float(vector_size,1)

wgrgr_boast = ANArray.float(64,vector_size,pnode,pnode,pgaus)
agrau_boast = ANArray.float(64,vector_size,pnode,pgaus)
elauu_boast = ANArray.float(64,vector_size,pnode*ndime,pnode*ndime)
elaup_boast = ANArray.float(64,vector_size,pnode*ndime,pnode)
elapp_boast = ANArray.float(64,vector_size,pnode,pnode)
elapu_boast = ANArray.float(64,vector_size,pnode,pnode*ndime)
elrbu_boast = ANArray.float(64,vector_size,ndime,pnode)
elrbp_boast = ANArray.float(64,vector_size,pnode)

elauq_boast = ANArray.float(vector_size,pnode*ndime,1)
elapq_boast = ANArray.float(vector_size,pnode,1)
elaqu_boast = ANArray.float(vector_size,1,pnode*ndime)
elaqp_boast = ANArray.float(vector_size,1,pnode)
elaqq_boast = ANArray.float(vector_size,1,1)
elrbq_boast = ANArray.float(vector_size,1)

stats_ref[i] = k_orig.run(pnode,pgaus,gpden,gpvis,gppor,gpsp1,gpsp2,gpvol,gpsha,gpcar,gpadv,
                           gpvep,gpgrp,gprhs,gprhc,gpvel,gpsgs,elvel,elpre,elbub,elauu_ref,elaup_ref,
                           elapp_ref,elapu_ref,elrbu_ref,elrbp_ref,dtinv_loc,dtsgs,pbubl,gpsha_bub,gpcar_bub,
                           elauq_ref,elapq_ref,elaqu_ref,elaqp_ref,elaqq_ref,elrbq_ref,
                           # Original global variables
                           kfl_lumped,mnode,ntens,kfl_duatss,fact_duatss,kfl_stabi_nsi,fvins_nsi,
                           fcons_nsi,bemol_nsi,kfl_regim_nsi,fvela_nsi,kfl_rmom2_nsi,kfl_press_nsi,
                           kfl_p1ve2_nsi,kfl_linea_nsi,kfl_confi_nsi,nbdfp_nsi,
                           kfl_sgsti_nsi,kfl_nota1_nsi,kfl_limit_nsi,kfl_penal_nsi,penal_nsi,
                           kfl_bubbl_nsi,ndime,agrau_ref,wgrgr_ref)

stats_boast[i] = k_boast.run(pnode,pgaus,gpden,gpvis,gppor,gpsp1,gpsp2,gpvol,gpsha,gpcar,gpadv,
                           gpvep,gpgrp,gprhs,gprhc,gpvel,gpsgs,elvel,elpre,elbub,elauu_boast,elaup_boast,
                           elapp_boast,elapu_boast,elrbu_boast,elrbp_boast,dtinv_loc,dtsgs,pbubl,gpsha_bub,gpcar_bub,
                           elauq_boast,elapq_boast,elaqu_boast,elaqp_boast,elaqq_boast,elrbq_boast,
                           # Original global variables
                           kfl_lumped,mnode,ntens,kfl_duatss,fact_duatss,kfl_stabi_nsi,fvins_nsi,
                           fcons_nsi,bemol_nsi,kfl_regim_nsi,fvela_nsi,kfl_rmom2_nsi,kfl_press_nsi,
                           kfl_p1ve2_nsi,kfl_linea_nsi,kfl_confi_nsi,nbdfp_nsi,
                           kfl_sgsti_nsi,kfl_nota1_nsi,kfl_limit_nsi,kfl_penal_nsi,penal_nsi,
                           kfl_bubbl_nsi,ndime,agrau_boast,wgrgr_boast)

epsilon = 10e-15

diff_agrau = (agrau_ref - agrau_boast).abs
diff_wgrgr = (wgrgr_ref - wgrgr_boast).abs
diff_elauu = (elauu_ref - elauu_boast).abs
diff_elrbu = (elrbu_ref - elrbu_boast).abs
diff_elapu = (elapu_ref - elapu_boast).abs
diff_elaup = (elaup_ref - elaup_boast).abs
diff_elapp = (elapp_ref - elapp_boast).abs
diff_elrbp = (elrbp_ref - elrbp_boast).abs

raise "Error: residue too big for agrau" if (diff_agrau > epsilon).to_a.flatten.include? 1
raise "Error: residue too big for wgrgr" if (diff_wgrgr > epsilon).to_a.flatten.include? 1
raise "Error: residue too big for elauu" if (diff_elauu > epsilon).to_a.flatten.include? 1
raise "Error: residue too big for elrbu" if (diff_elrbu > epsilon).to_a.flatten.include? 1
raise "Error: residue too big for elapu" if (diff_elapu > epsilon).to_a.flatten.include? 1
raise "Error: residue too big for elaup" if (diff_elaup > epsilon).to_a.flatten.include? 1
raise "Error: residue too big for elapp" if (diff_elapp > epsilon).to_a.flatten.include? 1
raise "Error: residue too big for elrbp" if (diff_elrbp > epsilon).to_a.flatten.include? 1

puts "Done"

# }
# t = []
# t[1] = stats_ref.inject(0){|sum,e| sum+e[:duration]} / stats_ref.length
# t[2] = stats_boast.inject(0){|sum,e| sum+e[:duration]} / stats_boast.length
#  t
# wgrgr_ref
