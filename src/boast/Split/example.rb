require 'BOAST'
include BOAST
require './KSplitOssRef_v1.rb'
require './KSplitOssRef_v2.rb'
require './KSplitOssBoast.rb'
require 'narray_ffi'
require '../Common/CommonArgs.rb'

class Example
  # To retrieve common arguments settings
  include CommonArgs    
  
  def self.run
    # Specify the kernel options
    nests = [1,2,3,4,5,6,7,8,9,10,11]
    vector_size=2
    ndime = 3
    k_boast_params = {:vector_length => vector_size, :nests => nests, :unroll => false, :inline => :inlined}
    
    set_lang(C)
    # Create the kernel, ndime is needed to create the BOAST version
    k = KSplitBoast::new(k_boast_params,ndime)
    # Generate the will create a CKernel instance available using the variable kernel.
    # And can be used directly used as usuall 
    k.generate
    puts k.kernel   
    k.kernel.build()
    
    # The CommonArgs module need to be initialized.
    # We need to specify the number of kernel instances 
    # this way the module can create the correct number 
    # of arguments.
    seed = 10
    pgaus = 8
    pnode = 8
    CommonArgs.init(pgaus,pnode,ndime,vector_size,1,seed)
    
    # Arguments can still be set manually especially for those which the path 
    # taken depends of their value
    @@kfl_lumped = 1
    @@kfl_limit_nsi = 1
    @@kfl_stabi_nsi = 1

    # Then just run the kernel using the reference to the arguments
    puts k.kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
                      @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
                      @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                      @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                      @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
                      @@gpsha,@@gpcar,@@gpadv,@@gpvep[0],@@gpgrp[0],@@gprhs[0],@@gprhc[0],@@gpvel,
                      @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[0],@@agrau[0],@@elauu[0],@@elaup[0],
                      @@elapp[0],@@elapu[0],@@elrbu[0],@@elrbp[0],@@dtinv_loc,@@dtsgs,@@pbubl,
                      @@gpsha_bub,@@gpcar_bub,@@elauq[0],@@elapq[0],@@elaqu[0],@@elaqp[0],
                      @@elaqq[0],@@elrbq[0])

  end
end
Example.run
