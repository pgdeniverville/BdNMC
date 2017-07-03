#include <vector>
#include <string>

using std::string;
using std::vector;

const string beam_energy_key = "beam_energy";
const string max_dm_energy_key = "max_dm_energy";
const string dm_energy_res_key = "dm_energy_resolution";
const string min_scatter_energy_key = "min_scatter_energy";
const string max_scatter_energy_key = "max_scatter_energy";
const string angle_lower_limit_key = "min_scatter_angle";
const string angle_upper_limit_key = "max_scatter_angle";
const string timing_key = "timing_cut";
const string scatter_dist_filename_key = "scatter_dist_filename";

const string production_distribution_key = "production_distribution";
const string particle_list_file_key = "particle_list_file";
const string prodkey = "production_channel";
const string meson_per_pi0_key = "meson_per_pi0";
const string parton_V_neutron_file_key = "parton_V_neutron_file";
const string parton_V_proton_file_key = "parton_V_proton_file";
const string sanford_wang_key = "sanfordwang_file";
const string distribution_parameter_key = "distribution_parameter_file";
const string ptmax_key = "ptmax";
const string ptmin_key = "ptmin";
const string zmin_key = "zmin";
const string zmax_key = "zmax";
const string part_list_pos_key = "particle_list_position";

const string dist_mod_key = "dist_mod";
const string pos_offset_key = "position_offset";

const string detkey = "detector";

const string modelkey = "model";

const string dmkey = "dark_matter_mass";
const string Vkey = "dark_photon_mass";
const string epskey = "epsilon";
const string alkey = "alpha_D";
const string signalkey = "signal_channel";
const string coherent_key = "coherent";
const string outkey = "output_file";
const string sumkey = "summary_file";
const string output_mode_key = "output_mode"; const string out_mode_def = "comprehensive";
const string efficiency_key = "efficiency";
const string POTkey = "POT";
const string pi0_ratio_key = "pi0_per_POT";
const string burn_in_key = "burn_max";
const string burn_timeout_key = "burn_timeout";
const string repeat_key = "repeat_count";
const string max_trials_key = "max_trials";
const string seed_key = "seed";
const string p_cross_key = "proton_target_cross_section";

vector<string> spherearr = {"radius","x-position","y-position","z-position"};
vector<string> cylinderarr = {"radius","length","x-position","y-position","z-position","det-theta","det-phi"};
vector<string> cuboidarr = {"length","width","height","x-position","y-position","z-position","det-phi","det-theta","det-psi"};
const string run_key = "run";
