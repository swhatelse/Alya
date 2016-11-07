module Params
  def self.init(vector_size,dim,seed)
    ANArray.srand(seed) if seed

    @@mnode = 10
    @@ntens = 2
    @@kfl_duatss = 1
    @@fact_duatss = 2
    @@fvins_nsi = 2.0
    @@fcons_nsi = 0.5
    @@bemol_nsi = 1.0
    @@kfl_regim_nsi = 3
    @@fvela_nsi = ANArray.float(64,3).random!
    @@kfl_rmom2_nsi = 1
    @@kfl_press_nsi = 1
    @@kfl_p1ve2_nsi = 1
    @@kfl_linea_nsi = 2

    @@kfl_confi_nsi = 1
    @@nbdfp_nsi = 5
    @@kfl_sgsti_nsi = 1
    @@kfl_nota1_nsi = 0
    @@kfl_penal_nsi = 1
    @@penal_nsi = 1.0
    @@kfl_bubbl_nsi = 1

    @@pnode = 10
    @@pgaus = 10
    @@pevat = 10
    @@ndime = 3
    @@fvins_nsi = 2.0
    @@kfl_limit_nsi = 1
    @@kfl_sgsti_nsi = 1
    @@nbdfp_nsi = 3

    @@gpden = ANArray.float(64,vector_size,@@pgaus).random!
    @@gpvis = ANArray.float(64,vector_size,@@pgaus).random!
    @@gppor = ANArray.float(64,vector_size,@@pgaus).random! 
    @@gpsp1 = ANArray.float(64,vector_size,@@pgaus).random! 
    @@gpsp2 = ANArray.float(64,vector_size,@@pgaus).random! 
    @@gpvol = ANArray.float(64,vector_size,@@pgaus).random! 
    @@gpsha = ANArray.float(64,vector_size,@@pnode,@@pgaus).random! 
    @@gpcar = ANArray.float(64,vector_size,@@ndime,@@mnode,@@pgaus).random!
    @@gpadv = ANArray.float(64,vector_size,@@ndime,@@pgaus).random! 
    @@gpprp = ANArray.float(64,vector_size,@@pgaus).random! 
    @@gpvel = ANArray.float(64,vector_size,@@ndime,@@pgaus,10).random!
    @@gpsgs = ANArray.float(64,vector_size,@@ndime,@@pgaus,10).random!
    @@elvel = ANArray.float(64,vector_size,@@ndime,@@pnode,10).random!
    @@elpre = ANArray.float(64,vector_size,@@pnode,10).random!
    @@elbub = ANArray.float(64,vector_size).random!

    @@dtinv_loc = ANArray.float(64,vector_size).random!
    @@dtsgs = ANArray.float(64,vector_size).random!
    @@pbubl = ANArray.int(64,vector_size).random!
    @@gpsha_bub = ANArray.float(64,vector_size,@@pgaus).random!
    @@gpcar_bub = ANArray.float(64,vector_size,@@ndime,@@pgaus).random!

    @@wgrgr_ref = ANArray.float(64,vector_size,@@pnode,@@pnode,@@pgaus)
    @@agrau_ref = ANArray.float(64,vector_size,@@pnode,@@pgaus)
    @@elauu_ref = ANArray.float(64,vector_size,@@pnode*@@ndime,@@pnode*@@ndime)
    @@elaup_ref = ANArray.float(64,vector_size,@@pnode*@@ndime,@@pnode)
    @@elapp_ref = ANArray.float(64,vector_size,@@pnode,@@pnode)
    @@elapu_ref = ANArray.float(64,vector_size,@@pnode,@@pnode*@@ndime)
    @@elrbu_ref = ANArray.float(64,vector_size,@@ndime,@@pnode)
    @@elrbp_ref = ANArray.float(64,vector_size,@@pnode)

    @@elauq_ref = ANArray.float(vector_size,@@pnode*@@ndime,1)
    @@elapq_ref = ANArray.float(vector_size,@@pnode,1)
    @@elaqu_ref = ANArray.float(vector_size,1,@@pnode*@@ndime)
    @@elaqp_ref = ANArray.float(vector_size,1,@@pnode)
    @@elaqq_ref = ANArray.float(vector_size,1,1)
    @@elrbq_ref = ANArray.float(vector_size,1)

    @@gpgrp_ref = ANArray.float(64,vector_size,@@ndime,@@pgaus).random! 
    @@gprhs_ref = ANArray.float(64,vector_size,@@ndime,@@pgaus).random!
    @@gpvep_ref = ANArray.float(64,vector_size,@@ndime,@@pgaus).random! 
    @@gprhc_ref = ANArray.float(64,vector_size,@@pgaus).random! 

    @@wgrgr_boast = ANArray.float(64,vector_size,@@pnode,@@pnode,@@pgaus)
    @@agrau_boast = ANArray.float(64,vector_size,@@pnode,@@pgaus)
    @@elauu_boast = ANArray.float(64,vector_size,@@pnode*@@ndime,@@pnode*@@ndime)
    @@elaup_boast = ANArray.float(64,vector_size,@@pnode*@@ndime,@@pnode)
    @@elapp_boast = ANArray.float(64,vector_size,@@pnode,@@pnode)
    @@elapu_boast = ANArray.float(64,vector_size,@@pnode,@@pnode*@@ndime)
    @@elrbu_boast = ANArray.float(64,vector_size,@@ndime,@@pnode)
    @@elrbp_boast = ANArray.float(64,vector_size,@@pnode)

    @@elauq_boast = ANArray.float(vector_size,@@pnode*@@ndime,1)
    @@elapq_boast = ANArray.float(vector_size,@@pnode,1)
    @@elaqu_boast = ANArray.float(vector_size,1,@@pnode*@@ndime)
    @@elaqp_boast = ANArray.float(vector_size,1,@@pnode)
    @@elaqq_boast = ANArray.float(vector_size,1,1)
    @@elrbq_boast = ANArray.float(vector_size,1)

    @@gpgrp_boast = @@gpgrp_ref.clone
    @@gprhs_boast = @@gprhs_ref.clone
    @@gpvep_boast = @@gpvep_ref.clone
    @@gprhc_boast = @@gprhc_ref.clone
  end
end
