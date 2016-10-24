gem "minitest"
require 'minitest/autorun'
require 'BOAST'
include BOAST
require '/tmp/alya.rb'

class TestFor < Minitest::Test
  def init(vector_size,seed)
    NArray.srand(seed) if seed
    @pnode = 2
    @mnode = 2
    @pgaus = 2
    @pevat = 2
    @ndime = 3
    @fvins_nsi = 1.0
    @kfl_lumped = 1
    @kfl_limit_nsi = 1
    @kfl_sgsti_nsi = 1
    @nbdfp_nsi = 3

    @gpden = NArray.float(vector_size,@pgaus).random
    @gpvis = NArray.float(vector_size,@pgaus).random
    @gppor = NArray.float(vector_size,@pgaus).random 
    @gpsp1 = NArray.float(vector_size,@pgaus).random 
    @gpsp2 = NArray.float(vector_size,@pgaus).random 
    @gpvol = NArray.float(vector_size,@pgaus).random 
    @gpsha = NArray.float(vector_size,@pnode,@pgaus).random 
    @gpcar = NArray.float(vector_size,@ndime,@mnode,@pgaus).random
    @gpadv = NArray.float(vector_size,@ndime,@pgaus).random 
    @gpvep = NArray.float(vector_size,@ndime,@pgaus).random 
    @gpprp = NArray.float(vector_size,@pgaus).random 
    @gpgrp = NArray.float(vector_size,@ndime,@pgaus).random 
    @gprhs = NArray.float(vector_size,@ndime,@pgaus).random 
    @gpvel = NArray.float(vector_size,@ndime,@pgaus,10).random
    @gpsgs = NArray.float(vector_size,@ndime,@pgaus,10).random
    @elvel = NArray.float(vector_size,@ndime,@pnode,10).random
    @dtinv_loc = NArray.float(vector_size).random
    @dtsgs = NArray.float(vector_size).random

    @wgrgr_ref = NArray.float(vector_size,@pnode,@pnode,@pgaus)
    @agrau_ref = NArray.float(vector_size,@pnode,@pgaus)
    @elauu_ref = NArray.float(vector_size,@pnode*@ndime,@pnode*@ndime)
    @elaup_ref = NArray.float(vector_size,@pnode*@ndime,@pnode)
    @elapp_ref = NArray.float(vector_size,@pnode,@pnode)
    @elapu_ref = NArray.float(vector_size,@pnode,@pnode*@ndime)
    @elrbu_ref = NArray.float(vector_size,@ndime,@pnode)
    @elrbp_ref = NArray.float(vector_size,@pnode)

    @wgrgr_boast = NArray.float(vector_size,@pnode,@pnode,@pgaus)
    @agrau_boast = NArray.float(vector_size,@pnode,@pgaus)
    @elauu_boast = NArray.float(vector_size,@pnode*@ndime,@pnode*@ndime)
    @elaup_boast = NArray.float(vector_size,@pnode*@ndime,@pnode)
    @elapp_boast = NArray.float(vector_size,@pnode,@pnode)
    @elapu_boast = NArray.float(vector_size,@pnode,@pnode*@ndime)
    @elrbu_boast = NArray.float(vector_size,@ndime,@pnode)
    @elrbp_boast = NArray.float(vector_size,@pnode)
  end
  def test_nests_orig_vs_boast
    vector_size=2

    k_orig_params = {:vector_length => vector_size, :preprocessor => false, :nests => [1]}
    k_boast_params = {:vector_length => vector_size, :nests => [1]}

    k_orig = generate_ref(k_orig_params)
    k_orig.build(:FCFLAGS => "-cpp")

    set_lang(FORTRAN)
    k_boast = generate_boast_implem(k_boast_params)
    k_boast.build

    init(vector_size,10)

    k_orig.run(@ndime,@mnode,@pnode,@pgaus,@pevat,@gpden,@gpvis,@gppor,
               @gpsp1,@gpsp2,@gpvol,@gpsha,@gpcar,@gpadv,@gpvep,@gpprp,
               @gpgrp,@gprhs,@gpvel,@gpsgs,@wgrgr_ref,@agrau_ref,@elvel,@elauu_ref,
               @elaup_ref,@elapp_ref,@elapu_ref,@elrbu_ref,@elrbp_ref,@dtinv_loc,@dtsgs,
               @fvins_nsi,@kfl_lumped,@kfl_limit_nsi,@kfl_sgsti_nsi,@nbdfp_nsi)

    k_boast.run(@ndime,@mnode,@pnode,@pgaus,@pevat,@gpden,@gpvis,@gppor,
                @gpsp1,@gpsp2,@gpvol,@gpsha,@gpcar,@gpadv,@gpvep,@gpprp,
                @gpgrp,@gprhs,@gpvel,@gpsgs,@wgrgr_boast,@agrau_boast,@elvel,@elauu_boast,
                @elaup_boast,@elapp_boast,@elapu_boast,@elrbu_boast,@elrbp_boast,@dtinv_loc,@dtsgs,
                @fvins_nsi,@kfl_lumped,@kfl_limit_nsi,@kfl_sgsti_nsi,@nbdfp_nsi)

    assert(@agrau_ref == @agrau_boast)
    assert(@wgrgr_ref == @wgrgr_boast)
  end

  def test_nests_fortran_vs_c
    vector_size=2

    k_orig_params = {:vector_length => vector_size, :preprocessor => false, :nests => [1]}
    k_boast_params = {:vector_length => vector_size, :nests => [1]}

    set_lang(FORTRAN)
    k_c = generate_boast_implem(k_boast_params)
    k_c.build

    set_lang(FORTRAN)
    k_fortran = generate_boast_implem(k_boast_params)
    k_fortran.build

    init(vector_size,10)

    k_c.run(@ndime,@mnode,@pnode,@pgaus,@pevat,@gpden,@gpvis,@gppor,
               @gpsp1,@gpsp2,@gpvol,@gpsha,@gpcar,@gpadv,@gpvep,@gpprp,
               @gpgrp,@gprhs,@gpvel,@gpsgs,@wgrgr_ref,@agrau_ref,@elvel,@elauu_ref,
               @elaup_ref,@elapp_ref,@elapu_ref,@elrbu_ref,@elrbp_ref,@dtinv_loc,@dtsgs,
               @fvins_nsi,@kfl_lumped,@kfl_limit_nsi,@kfl_sgsti_nsi,@nbdfp_nsi)

    k_fortran.run(@ndime,@mnode,@pnode,@pgaus,@pevat,@gpden,@gpvis,@gppor,
                @gpsp1,@gpsp2,@gpvol,@gpsha,@gpcar,@gpadv,@gpvep,@gpprp,
                @gpgrp,@gprhs,@gpvel,@gpsgs,@wgrgr_boast,@agrau_boast,@elvel,@elauu_boast,
                @elaup_boast,@elapp_boast,@elapu_boast,@elrbu_boast,@elrbp_boast,@dtinv_loc,@dtsgs,
                @fvins_nsi,@kfl_lumped,@kfl_limit_nsi,@kfl_sgsti_nsi,@nbdfp_nsi)

    assert(@agrau_ref == @agrau_boast)
    assert(@wgrgr_ref == @wgrgr_boast)
  end
end
