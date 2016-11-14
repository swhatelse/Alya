require '/tmp/alya.rb'
require 'narray_ffi'
require_relative './mod_set_params.rb'
require 'yaml'
require 'pp'
require 'csv'
require 'optparse'
require_relative './LogInfo.rb'

class Experiment
  include Params
  def self.run(output_info)
    # nests = (1..11).to_a
    LogInfo.init(output_info)
    LogInfo.get_info
    nests = [1,2,3,4,5,6,7,9]
    vector_size=2
    dimension=2
    seed = 10
    epsilon = 10e-15

    stats = {}
    k = {}
    nests.each{|n|
      nest_name = ":nest#{n}"
      opts = {:vector_length => vector_size, :preprocessor => false, :nests => [n], :unroll => true}
      h_ref = {:kernel => :ref, :nest => n}
      h_boast = {:kernel => :boast, :nest => n}

      set_fortran_line_length(100)
      k[h_ref] = generate_ref_v2(opts)
      LogInfo.register_kernel_info(h_ref, k[h_ref].to_s)
      k[h_ref].build(:FCFLAGS => "-O3")

      set_lang(C)
      k[h_boast] = generate_boast_implem(opts)
      LogInfo.register_kernel_info(h_boast, k[h_boast].to_s)
      k[h_boast].build(:CFLAGS => "-O3")
      
      stats[h_ref] = []
      stats[h_boast] = []

      Params.init(vector_size,dimension,seed)
      @@kfl_lumped = 2 # 1
      @@kfl_limit_nsi = 1 # 2
      @@kfl_stabi_nsi = 1 # -1
      100.times{|i|
        k[h_ref].run(@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,@@gpsha,@@gpcar,@@gpadv,
                   @@gpvep_ref,@@gpgrp_ref,@@gprhs_ref,@@gprhc_ref,@@gpvel,@@gpsgs,@@elvel,@@elpre,@@elbub,@@elauu_ref,@@elaup_ref,
                   @@elapp_ref,@@elapu_ref,@@elrbu_ref,@@elrbp_ref,@@dtinv_loc,@@dtsgs,@@pbubl,@@gpsha_bub,@@gpcar_bub,
                   @@elauq_ref,@@elapq_ref,@@elaqu_ref,@@elaqp_ref,@@elaqq_ref,@@elrbq_ref,
                   # Original global variables
                   @@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,@@fvins_nsi,
                   @@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,@@kfl_press_nsi,
                   @@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                   @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                   @@kfl_bubbl_nsi,@@ndime,@@agrau_ref,@@wgrgr_ref)

        k[h_boast].run(@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,@@gpsha,@@gpcar,@@gpadv,
                    @@gpvep_boast,@@gpgrp_boast,@@gprhs_boast,@@gprhc_boast,@@gpvel,@@gpsgs,@@elvel,@@elpre,@@elbub,@@elauu_boast,@@elaup_boast,
                    @@elapp_boast,@@elapu_boast,@@elrbu_boast,@@elrbp_boast,@@dtinv_loc,@@dtsgs,@@pbubl,@@gpsha_bub,@@gpcar_bub,
                    @@elauq_boast,@@elapq_boast,@@elaqu_boast,@@elaqp_boast,@@elaqq_boast,@@elrbq_boast,
                    # Original global variables
                    @@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,@@fvins_nsi,
                    @@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,@@kfl_press_nsi,
                    @@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                    @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                    @@kfl_bubbl_nsi,@@ndime,@@agrau_boast,@@wgrgr_boast)
      }
      
      1000.times{|i|
        Params.init(vector_size,dimension,seed)
        @@kfl_lumped = 2 # 1
        @@kfl_limit_nsi = 1 # 2
        @@kfl_stabi_nsi = 1 # -1
        
        stats[h_ref][i] = k[h_ref].run(@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,@@gpsha,@@gpcar,@@gpadv,
                                                  @@gpvep_ref,@@gpgrp_ref,@@gprhs_ref,@@gprhc_ref,@@gpvel,@@gpsgs,@@elvel,@@elpre,@@elbub,@@elauu_ref,@@elaup_ref,
                                                  @@elapp_ref,@@elapu_ref,@@elrbu_ref,@@elrbp_ref,@@dtinv_loc,@@dtsgs,@@pbubl,@@gpsha_bub,@@gpcar_bub,
                                                  @@elauq_ref,@@elapq_ref,@@elaqu_ref,@@elaqp_ref,@@elaqq_ref,@@elrbq_ref,
                                                  # Original global variables
                                                  @@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,@@fvins_nsi,
                                                  @@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,@@kfl_press_nsi,
                                                  @@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                                                  @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                                                  @@kfl_bubbl_nsi,@@ndime,@@agrau_ref,@@wgrgr_ref)[:duration]

        stats[h_boast][i] = k[h_boast].run(@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,@@gpsha,@@gpcar,@@gpadv,
                                                     @@gpvep_boast,@@gpgrp_boast,@@gprhs_boast,@@gprhc_boast,@@gpvel,@@gpsgs,@@elvel,@@elpre,@@elbub,@@elauu_boast,@@elaup_boast,
                                                     @@elapp_boast,@@elapu_boast,@@elrbu_boast,@@elrbp_boast,@@dtinv_loc,@@dtsgs,@@pbubl,@@gpsha_bub,@@gpcar_bub,
                                                     @@elauq_boast,@@elapq_boast,@@elaqu_boast,@@elaqp_boast,@@elaqq_boast,@@elrbq_boast,
                                                     # Original global variables
                                                     @@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,@@fvins_nsi,
                                                     @@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,@@kfl_press_nsi,
                                                     @@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                                                     @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                                                     @@kfl_bubbl_nsi,@@ndime,@@agrau_boast,@@wgrgr_boast)[:duration]

        diff_agrau = (@@agrau_ref - @@agrau_boast).abs
        diff_wgrgr = (@@wgrgr_ref - @@wgrgr_boast).abs
        diff_elauu = (@@elauu_ref - @@elauu_boast).abs
        diff_elrbu = (@@elrbu_ref - @@elrbu_boast).abs
        diff_elapu = (@@elapu_ref - @@elapu_boast).abs
        diff_elaqu = (@@elaqu_ref - @@elaqu_boast).abs
        diff_elaqp = (@@elaqp_ref - @@elaqp_boast).abs
        diff_elaqq = (@@elaqq_ref - @@elaqq_boast).abs
        diff_elapq = (@@elapq_ref - @@elapq_boast).abs
        diff_elauq = (@@elauq_ref - @@elauq_boast).abs
        diff_elaup = (@@elaup_ref - @@elaup_boast).abs
        diff_elapp = (@@elapp_ref - @@elapp_boast).abs
        diff_elrbp = (@@elrbp_ref - @@elrbp_boast).abs
        diff_elrbq = (@@elrbq_ref - @@elrbq_boast).abs

        diff_gpgrp = (@@gpgrp_boast - @@gpgrp_ref).abs
        diff_gprhs = (@@gprhs_boast - @@gprhs_ref).abs
        diff_gpvep = (@@gpvep_boast - @@gpvep_ref).abs
        diff_gprhc = (@@gprhc_boast - @@gprhc_ref).abs


        raise "Error: residue too big for agrau" if (diff_agrau > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for wgrgr" if (diff_wgrgr > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elauu" if (diff_elauu > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elrbu" if (diff_elrbu > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elapu" if (diff_elapu > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elaqu" if (diff_elaqu > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elaqp" if (diff_elaqp > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elaqq" if (diff_elaqq > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elapq" if (diff_elapq > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elauq" if (diff_elauq > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elaup" if (diff_elaup > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elapp" if (diff_elapp > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elrbp" if (diff_elrbp > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elrbq" if (diff_elrbq > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for gpgrp" if (diff_gpgrp > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for gprhs" if (diff_gprhs > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for gprhc" if (diff_gprhc > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for gpvep" if (diff_gpvep > epsilon).to_a.flatten.include? 1
      }
    }
    LogInfo.dump_info
    puts "Done"
    return stats,k
  end
end

options = {}

opt_parser = OptionParser.new { |opts|
  opts.banner = "Usage: run.rb --[options]=[value]"

  opts.on("-dVAL", "--data=VAL", "Specify the path where to store the data") { |n|
    options[:data_path] = n
  }

  opts.on("-iVAL", "--info=VAL", "Specify the path where to store the information such as the generated source, compiler info, etc...") { |n|
    options[:info_path] = n
  }

  opts.on("-h", "--help", "Prints this help") {
    puts opts
    exit
  }
}.parse!

stats, infos = Experiment.run(options[:info_path])

if options[:data_path] then
  File::open( options[:data_path], "w") { |f|
    f.print YAML::dump(stats)
  }
end
