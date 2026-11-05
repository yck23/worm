#include "worm.hpp"

#define WORM_STRINGIFY_IMPL(value) #value
#define WORM_STRINGIFY(value) WORM_STRINGIFY_IMPL(value)


void worm::print_params(std::ostream& os) const {
#ifdef UNISYS
  os << "# config : UNIFORM                                  : yes\n";
#else
  os << "# config : UNIFORM                                  : no\n";
#endif
#ifdef CAN_WINDOW
  os << "# config : CWINDOW                                  : yes\n";
  os << "# canonical open-worm window                       : " << can_window << "\n";
#else
  os << "# config : CWINDOW                                  : no\n";
#endif
  Diagram_type::print_name(os);
  os << "# run time limit                                    : " << runtimelimit << "\n";
  os << "# statistics were force-reset on restore             : " << reset_statistics << "\n";
  os << "# maximum number of sweeps                          : " << total_sweeps << "\n";
  os << "# maximum number of sweeps for thermalization       : " << thermalization_sweeps << "\n";
  os << "# inverse temperature                               : " << beta << "\n";
  os << "# canonical measurement value                       : " << canonical << "\n";
  os << "# random number seed                                : " << parameters["seed"].as<size_t>() << "\n";
  os << "# Ntest                                             : " << Ntest<< "\n";
  os << "# Nsave                                             : " << Nsave<< "\n";
  os << "# Nmeasure O(1)                                     : " << Nmeasure << "\n";
  os << "# Nmeasure O(N)                                     : " << Nmeasure2 << "\n";
  os << "# Nresponse_measure (0 disables)                    : " << Nresponse_measure << "\n";
  os << "# Nresponse_max_blocks                              : " << Nresponse_max_blocks << "\n";
  os << "# energy shift                                      : " << E_off << "\n";
  os << "# weight worm configurations                        : " << C_worm << "\n";
  os << "# update probability insert worm                    : " << update_prob[update_tag::insertworm]    << "\t" << update_prob_cuml[update_tag::insertworm]    << "\n";
  os << "# update probability move worm                      : " << update_prob[update_tag::moveworm]      << "\t" << update_prob_cuml[update_tag::moveworm]    << "\n";
  os << "# update probability insert kink                    : " << update_prob[update_tag::insertkink]    << "\t" << update_prob_cuml[update_tag::insertkink]    << "\n";
  os << "# update probability delete kink                    : " << update_prob[update_tag::deletekink]    << "\t" << update_prob_cuml[update_tag::deletekink]    << "\n";
  os << "# update probability glue worm                      : " << update_prob[update_tag::glueworm]      << "\t" << update_prob_cuml[update_tag::glueworm]      << "\n";
  
  MyLatt->print_name(os);
  MyLatt->print_params(os);
  MyModel->print_params(os);
  os << "# on-site V : enabled                               : " << (onsite_potential_enabled ? "yes" : "no (V_i=0; configured file ignored)") << "\n";
  if (ModelClassifierName == "BoseHubbard" && !onsite_potential.empty()) {
    const auto bounds = std::minmax_element(onsite_potential.begin(), onsite_potential.end());
    const double mean = std::accumulate(onsite_potential.begin(), onsite_potential.end(), 0.) / onsite_potential.size();
    os << "# on-site V : source (provenance only on restore)     : " << onsite_potential_source << "\n";
    os << "# on-site V : min, mean, max                         : " << *bounds.first << "\t" << mean << "\t" << *bounds.second << "\n";
    os << "# on-site V : energy convention                     : +V_i*n_i - mu*n_i\n";
  }
  
}


boost::multi_array<double,2> worm::serialize_op_string() const {
  size_t maxsize = 0;
  for (auto const& vertexlist : operator_string){
    if (vertexlist.size() > maxsize)
      maxsize = vertexlist.size();
  }
  boost::multi_array<double, 2> op_string_arr(
    boost::extents[operator_string.size()][maxsize*(zcmax + 5) + 1]);

  auto row_it = op_string_arr.begin();

  for (auto const& vertexlist : operator_string){
    auto col_it = (row_it++)->begin();
    *(col_it++) = vertexlist.size();
    for (auto const& vertex : vertexlist){
      *(col_it++) = vertex.before();
      *(col_it++) = vertex.after();
      *(col_it++) = vertex.link();
      *(col_it++) = vertex.time();
      *(col_it++) = vertex.color();
      if (vertex.link() != -1) {
        for (size_t d = 0; d < zcmax; ++d, ++col_it)
          *col_it = vertex.get_assoc_time(d);
      }
    }
  }
  return op_string_arr;
}


void worm::save_hamiltonian(alps::hdf5::archive& ar, const std::string& path) const {
  const std::vector<int> shape = {
    static_cast<int>(MyLatt->get_Ls(0)), static_cast<int>(MyLatt->get_Ls(1)),
    static_cast<int>(MyLatt->get_Ls(2)), static_cast<int>(LATTICE::n_basis),
    static_cast<int>(latt_dim), static_cast<int>(Nsites)};
  const std::vector<int> periodic = {parameters["pbcx"].as<bool>(),
    parameters["pbcy"].as<bool>(), parameters["pbcz"].as<bool>()};
#ifdef UNISYS
  const int uniform_build = 1;
#else
  const int uniform_build = 0;
#endif
  ar[path + "/schema_version"] << 1;
  ar[path + "/model"] << ModelClassifierName;
  ar[path + "/lattice"] << std::string(WORM_STRINGIFY(LATTICE));
  ar[path + "/shape"] << shape;
  ar[path + "/shape_convention"] << std::string("Lx,Ly,Lz,n_basis,dimension,n_sites");
  ar[path + "/periodic"] << periodic;
  ar[path + "/uniform_build"] << uniform_build;
  ar[path + "/nmax"] << nmax;
  ar[path + "/onsite_potential"] << onsite_potential;
  ar[path + "/onsite_potential_enabled"] << static_cast<int>(onsite_potential_enabled);
  ar[path + "/source_file"] << onsite_potential_source;
  ar[path + "/site_index_convention"] << std::string("basis + n_basis * (x + Lx * (y + Ly * z))");
  ar[path + "/energy_convention"] << std::string("BoseHubbard: U*n*(n-1)/2 + onsite_potential*n - mu*n; Total_Energy includes -mu*N");
  ar[path + "/energy_units"] << std::string("same as the scalar Hamiltonian parameters; no implicit Hz conversion");
  const std::vector<std::string> names = ModelClassifierName == "BoseHubbard"
    ? std::vector<std::string>{"t_hop", "U_on", "mu", "V_nn"}
    : std::vector<std::string>{"Jpm", "Jzz", "h"};
  for (const auto& name : names) ar[path + "/" + name] << parameters[name].as<double>();
}

void worm::load_hamiltonian(alps::hdf5::archive& ar) {
  const std::string path = "checkpoint/hamiltonian";
  if (!ar.is_data(path + "/schema_version")) {
#ifndef UNISYS
    throw std::runtime_error("This checkpoint predates resolved nonuniform Hamiltonian storage; start a fresh programmable-potential run.");
#else
    return; // Historical uniform checkpoints already contain scalar parameters.
#endif
  }
  int version, stored_nmax, stored_uniform;
  std::string stored_model, stored_lattice;
  std::vector<int> stored_shape, stored_periodic;
  ar[path + "/schema_version"] >> version;
  ar[path + "/model"] >> stored_model;
  ar[path + "/lattice"] >> stored_lattice;
  ar[path + "/shape"] >> stored_shape;
  ar[path + "/periodic"] >> stored_periodic;
  ar[path + "/uniform_build"] >> stored_uniform;
  ar[path + "/nmax"] >> stored_nmax;
#ifdef UNISYS
  const int uniform_build = 1;
#else
  const int uniform_build = 0;
#endif
  const std::vector<int> shape = {
    static_cast<int>(MyLatt->get_Ls(0)), static_cast<int>(MyLatt->get_Ls(1)),
    static_cast<int>(MyLatt->get_Ls(2)), static_cast<int>(LATTICE::n_basis),
    static_cast<int>(latt_dim), static_cast<int>(Nsites)};
  const std::vector<int> periodic = {parameters["pbcx"].as<bool>(),
    parameters["pbcy"].as<bool>(), parameters["pbcz"].as<bool>()};
  if (version != 1 || stored_model != ModelClassifierName ||
      stored_lattice != WORM_STRINGIFY(LATTICE) || stored_shape != shape ||
      stored_periodic != periodic || stored_nmax != nmax || stored_uniform != uniform_build)
    throw std::runtime_error("Checkpoint Hamiltonian/lattice/build metadata does not match this simulation.");
  const std::vector<std::string> names = ModelClassifierName == "BoseHubbard"
    ? std::vector<std::string>{"t_hop", "U_on", "mu", "V_nn"}
    : std::vector<std::string>{"Jpm", "Jzz", "h"};
  for (const auto& name : names) {
    double stored;
    ar[path + "/" + name] >> stored;
    if (!std::isfinite(stored) || stored != parameters[name].as<double>())
      throw std::runtime_error("Checkpoint Hamiltonian parameter mismatch: " + name);
  }
  // This field is optional for schema-1 checkpoints written before the
  // explicit switch existed. The constructor already applies the old default.
  int stored_enabled = static_cast<int>(onsite_potential_enabled);
  if (ar.is_data(path + "/onsite_potential_enabled"))
    ar[path + "/onsite_potential_enabled"] >> stored_enabled;
  if (stored_enabled != static_cast<int>(onsite_potential_enabled))
    throw std::runtime_error("Checkpoint on-site potential enabled flag does not match the parameters.");
  ar[path + "/onsite_potential"] >> onsite_potential;
  ar[path + "/source_file"] >> onsite_potential_source;
  if (onsite_potential.size() != static_cast<std::size_t>(Nsites))
    throw std::runtime_error("Checkpoint on-site potential has the wrong number of sites.");
  for (double value : onsite_potential) {
    if (!std::isfinite(value)) throw std::runtime_error("Checkpoint on-site potential contains a non-finite value.");
    if (!onsite_potential_enabled && value != 0.)
      throw std::runtime_error("Checkpoint has a nonzero on-site potential while onsite_potential_enabled is false.");
  }
  MyModel->set_onsite_potential(onsite_potential);
}

void worm::save(alps::hdf5::archive & ar) const {
  // Most of the save logic is already implemented in the base class
  alps::mcbase::save(ar);
  save_hamiltonian(ar, "checkpoint/hamiltonian");

  // random number engine
  std::ostringstream engine_ss;
  engine_ss << MyGenerator;

  std::vector<double> wind_vec(mWinding.data(), mWinding.data() + mWinding.size());

  ar["checkpoint/random"] << engine_ss.str();
  ar["checkpoint/sweeps"] << sweeps;
  ar["checkpoint/counter"] << counter;
  ar["checkpoint/update_statistics"] << update_statistics;
  ar["checkpoint/configuration/operator_string"] << serialize_op_string();
  ar["checkpoint/configuration/mWinding"] << wind_vec;
  ar["checkpoint/configuration/nrvertex"] << nrvertex;
  ar["checkpoint/configuration/nr_of_particles"] << number_of_particles;
  ar["checkpoint/configuration/state"] << state;
  if (Nresponse_measure > 0) {
    ar["checkpoint/response/schema_version"] << 1;
    ar["checkpoint/response/interval"] << Nresponse_measure;
    ar["checkpoint/response/max_blocks"] << Nresponse_max_blocks;
    ar["checkpoint/response/partial_count"] << response_partial_count;
    ar["checkpoint/response/offsite_sum"] << response_offsite_sum;
    ar["checkpoint/response/density_sum"] << response_density_sum;
    ar["checkpoint/response/w2sum"] << response_w2sum;
    ar["checkpoint/response/chi_s"] << response_chi_s;
    ar["checkpoint/response/w2"] << response_w2;
    ar["checkpoint/response/density"] << response_density;
    ar["checkpoint/response/end_sweep"] << response_end_sweep;
    ar["checkpoint/response/sweep_count"] << response_sweep_count;
  }
  ar["checkpoint/configuration/hist_susceptibility_offsite"] << hist_susceptibility_offsite;
#ifdef UNISYS
  ar["checkpoint/configuration/hist_densmat"] << hist_densmat;
#ifdef MATSUBARA_MEAS
  if (Nfreq > 0) {
    ar["checkpoint/configuration/hist_grtau_re"] << hist_grtau_re;
    ar["checkpoint/configuration/hist_grtau_im"] << hist_grtau_im;
    // k-point and frequency metadata (deterministic, but stored for convenience)
    const size_t Lx = MyLatt->get_Ls(0);
    const size_t Ly = (LATTICE::dim > 1 ? MyLatt->get_Ls(1) : 1);
    vector<double> kpoints(2 * (size_t)Nsites);
    for (size_t r = 0; r < (size_t)Nsites; r++) {
      kpoints[2*r]   = 2.0 * M_PI * (r % Lx) / Lx;
      kpoints[2*r+1] = 2.0 * M_PI * (r / Lx) / Ly;
    }
    vector<double> freqs(Nfreq);
    for (size_t n = 0; n < Nfreq; n++) freqs[n] = 2.0 * M_PI * n / beta;
    ar["checkpoint/greenfun_kw/kpoints"] << kpoints;
    ar["checkpoint/greenfun_kw/matsubara_freqs"] << freqs;
    ar["checkpoint/greenfun_kw/flat_index_convention"] << std::string("k_idx * Nfreq + n");
  }
#endif
#endif

}


void worm::load(alps::hdf5::archive & ar) {
  // Validate against the parameters used to construct MyModel before ALPS
  // replaces parameters with the rank-local checkpoint's copy.
  load_hamiltonian(ar);
  alps::mcbase::load(ar);

  worm_diag = true;
  worm_at_stop = 0;
  worm_passes_nb_kink = 0;
  worm_dtime = 0.;
  new_measurement = true;
  worm_meas_densmat = false;


  auto find_elem = [&](Diagram_type::iterator begin, Diagram_type::iterator end, 
                       double findTime) -> Diagram_type::iterator
  {
    for(auto it = begin; it != end; ++it) {
        //if (it->time() == findTime)
        //    return it;
        if (std::abs(it->time() - findTime) <dtol) return it; 
    }
    return end;
  };

  std::string engine_str;
  std::vector<double> wind_vec;
  boost::multi_array<double, 2> op_string_arr;

  ar["checkpoint/random"] >> engine_str;
  ar["checkpoint/sweeps"] >> sweeps;
  ar["checkpoint/counter"] >> counter;
  ar["checkpoint/update_statistics"] >> update_statistics;
  ar["checkpoint/configuration/operator_string"] >> op_string_arr;
  ar["checkpoint/configuration/mWinding"] >> wind_vec;
  ar["checkpoint/configuration/nrvertex"] >> nrvertex;
  ar["checkpoint/configuration/nr_of_particles"] >> number_of_particles;
  ar["checkpoint/configuration/state"] >> state;

  // Old checkpoints have neither this group nor these parameter keys. Their
  // default remains disabled; an explicit parameter migration can enable a
  // fresh response record without borrowing any old density-matrix hits.
  const size_t restored_interval = parameters.exists("Nresponse_measure")
      ? parameters["Nresponse_measure"].as<size_t>() : 0;
  const size_t restored_limit = parameters.exists("Nresponse_max_blocks")
      ? parameters["Nresponse_max_blocks"].as<size_t>() : 1000000;
  if (restored_interval != Nresponse_measure || restored_limit != Nresponse_max_blocks)
    throw std::runtime_error("Response parameters changed between construction and checkpoint load.");
  clear_response_history();
  initialize_response_measurements();
  if (ar.is_group("checkpoint/response")) {
    if (Nresponse_measure == 0)
      throw std::runtime_error("Checkpoint contains response history but Nresponse_measure is disabled; refusing to discard it silently.");
    int schema_version = 0;
    size_t stored_interval = 0;
    ar["checkpoint/response/schema_version"] >> schema_version;
    ar["checkpoint/response/interval"] >> stored_interval;
    if (schema_version != 1 || stored_interval != Nresponse_measure)
      throw std::runtime_error("Unsupported response schema or changed response interval; use an explicit statistics migration.");
    ar["checkpoint/response/partial_count"] >> response_partial_count;
    ar["checkpoint/response/offsite_sum"] >> response_offsite_sum;
    ar["checkpoint/response/density_sum"] >> response_density_sum;
    ar["checkpoint/response/w2sum"] >> response_w2sum;
    ar["checkpoint/response/chi_s"] >> response_chi_s;
    ar["checkpoint/response/w2"] >> response_w2;
    ar["checkpoint/response/density"] >> response_density;
    ar["checkpoint/response/end_sweep"] >> response_end_sweep;
    ar["checkpoint/response/sweep_count"] >> response_sweep_count;
    const size_t blocks = response_chi_s.size();
    if (blocks > Nresponse_max_blocks || response_w2.size() != blocks ||
        response_density.size() != blocks || response_end_sweep.size() != blocks ||
        response_sweep_count.size() != blocks || response_partial_count >= Nresponse_measure ||
        !std::isfinite(response_offsite_sum) || response_offsite_sum < 0. ||
        !std::isfinite(response_density_sum) || response_density_sum < 0. ||
        !std::isfinite(response_w2sum) || response_w2sum < 0.)
      throw std::runtime_error("Invalid response checkpoint history or partial block.");
    for (size_t b = 0; b < blocks; ++b) {
      if (!std::isfinite(response_chi_s[b]) || !std::isfinite(response_w2[b]) ||
          !std::isfinite(response_density[b]) || response_chi_s[b] < 0. ||
          response_w2[b] < 0. || response_density[b] < 0. ||
          response_sweep_count[b] != Nresponse_measure || response_end_sweep[b] > sweeps ||
          (b > 0 && response_end_sweep[b] <= response_end_sweep[b - 1]))
        throw std::runtime_error("Invalid ordered response block in checkpoint.");
    }
    if (measurements["Superfluid_susceptibility"].count() != blocks ||
        measurements["Winding_number_squared_sum"].count() != blocks)
      throw std::runtime_error("Response accumulator counts disagree with the saved block history.");
    if (response_partial_count == 0 &&
        (response_offsite_sum != 0. || response_density_sum != 0. || response_w2sum != 0.))
      throw std::runtime_error("Empty response partial block contains nonzero sums.");
  }
  else if (Nresponse_measure > 0) {
    if (measurements["Superfluid_susceptibility"].count() != 0 ||
        measurements["Winding_number_squared_sum"].count() != 0)
      throw std::runtime_error("Response accumulators exist without their checkpoint history.");
    std::cout << "# No response history in legacy checkpoint; starting a new empty response record.\n";
  }
  else if ((measurements.has("Superfluid_susceptibility") &&
            measurements["Superfluid_susceptibility"].count() != 0) ||
           (measurements.has("Winding_number_squared_sum") &&
            measurements["Winding_number_squared_sum"].count() != 0)) {
    throw std::runtime_error("Disabled response measurements have stored samples without history.");
  }
  if (ar.is_data("checkpoint/configuration/hist_susceptibility_offsite"))
    ar["checkpoint/configuration/hist_susceptibility_offsite"] >> hist_susceptibility_offsite;
  else {
    hist_susceptibility_offsite = 0.;
    // Older disorder checkpoints have no equal-time histogram to restore.
    // Start a complete new measurement interval rather than retaining a partial one.
#ifndef UNISYS
    counter[counter_tag::MEASURE2] = 0;
#endif
  }
#ifdef UNISYS
  ar["checkpoint/configuration/hist_densmat"] >> hist_densmat;
  if (!ar.is_data("checkpoint/configuration/hist_susceptibility_offsite"))
    for (size_t i = 1; i < hist_densmat.size(); ++i)
      hist_susceptibility_offsite += hist_densmat[i];
#ifdef MATSUBARA_MEAS
  if (Nfreq > 0) {
    ar["checkpoint/configuration/hist_grtau_re"] >> hist_grtau_re;
    ar["checkpoint/configuration/hist_grtau_im"] >> hist_grtau_im;
  }
#endif
#endif

  std::istringstream engine_ss(engine_str);
  engine_ss >> MyGenerator;

  operator_string.clear();
  for (auto const& row : op_string_arr){
    Diagram_type vertexlist;
    auto col_it = row.begin() + 1;
    auto end = col_it + *row.begin() * (zcmax + 5);
    while (col_it != end) {
      auto before = StateType(*(col_it++));
      auto after = StateType(*(col_it++));
      auto link = SiteIndex(*(col_it++));
      auto time = *(col_it++);
      auto color = int(*(col_it++));
      vertexlist.insert(vertexlist.end(), Element(before, after, link, time, color));
      col_it += zcmax;
    }
    operator_string.push_back(vertexlist);
  }
  Nprtcls = 0;
  for (SiteIndex i = 0; i < Nsites; i++) {
    dummy_it[i] = find_elem(operator_string[i].begin(), operator_string[i].end(), beta);
    Nprtcls += dummy_it[i]->before() * beta;
  }
  dummy_it[Nsites] = find_elem(operator_string[Nsites].begin(), operator_string[Nsites].end(), beta);

  std::cout << "# Restoring existing associations ... ";
  for (SiteIndex i = 0; i < Nsites; i++) {
    for (size_t j = 0; j < zcmax; j++) {
      if (MyLatt->has_neighbor(i,j))
        dummy_it[i]->set_assoc(j, dummy_it[nb[i][j]]);
      else
        dummy_it[i]->set_assoc(j, *(dummy_it.end() - 1));
    }
  }

  for (SiteIndex i = 0; i < Nsites; i++) {
    size_t elem_idx = 0;
    for (auto elem_it = operator_string[i].begin(); elem_it != operator_string[i].end(); ++elem_it, ++elem_idx) {

      std::vector<double> assoc_vec(zcmax);
      auto row_it = op_string_arr.begin() + i;
      auto col_it = row_it->begin() + 1 + elem_idx*(zcmax + 5) + 5;
      std::copy(col_it, col_it + zcmax, assoc_vec.begin());

      for (size_t j = 0; j < zcmax; j++) {
        if (MyLatt->has_neighbor(i,j)) {
            auto assoc_elem_it = find_elem(operator_string[nb[i][j]].begin(), operator_string[nb[i][j]].end(), 
                                           assoc_vec.at(j));
            if (assoc_elem_it != operator_string[nb[i][j]].end())
                elem_it->set_assoc(j, assoc_elem_it);
            else
                throw std::runtime_error("Element not found during reconstruction of associations.");
        }
        else {
            elem_it->set_assoc(j, *(dummy_it.end() - 1));
        }
      }
    }
  }

  mWinding = Eigen::Map<LATTICE::Direction>(wind_vec.data());

  worm_head_it = dummy_it[0];
  worm_tail_it = dummy_it[0];

  std::cout << "...done. Computing potential energies...";
  Epot_tot = calc_potential_energy_nb() + calc_potential_energy_loc();
  Epot_measure = calc_potential_energy_measure();
  std::cout << "...done.\n";
  cout << "# Potential Energy tot : " << Epot_tot << endl;
  std::cout << "\n# Finished loading.\n";
  // Statistics are NOT reset here: an alps::params value read at this
  // point (post-restore) can never reflect anything other than what was
  // already stored in the checkpoint itself, since define_parameters()
  // (and therefore any CLI/INI override) is skipped whenever
  // parameters.is_restored() is true -- this is exactly the situation
  // load() is called in. See docs/reset_statistics_on_restore.md: the
  // caller (worm.run.cpp / worm.run_mpi.cpp) checks a separate,
  // manually-parsed --reset-statistics CLI flag and calls
  // force_reset_statistics() explicitly after load() returns, if given.
}

void worm::print_conf(std::ostream& os) const {
  os << "\n\n Printing operator string" << std::setprecision(16);
  for (SiteType i = 0; i < Nsites; i++) {
    os << "\nSite : " << i;
    int n = operator_string[i].size();
    int ii = 0;
    for (Diagram_type::const_iterator it = operator_string[i].begin(); it != operator_string[i].end(); ++it, ++ii) {
      it->print();
      if (ii >= n) {
        os << "\n Error with list!\n";
        throw std::runtime_error("print_conf: operator_string iteration exceeded list size");
      }
    }
    os << "\n--------------------\n\n\n";
  }
  os << "\n Worm head iterator time : "<< worm_head_it->time()  << "\n";
  os << "\n Worm tail iterator time : "<< worm_tail_it->time()  << "\n";
  os << "\n worm_at_stop = " << worm_at_stop << "\t worm_passes_nb_kink = " << worm_passes_nb_kink << "\t worm_measdensmat = " << worm_meas_densmat << "\n";
}


void worm::print_update_statistics(std::ostream& os) const {
  std::vector<string> name;
  os << setprecision(10);
  name.push_back("INSERT WORM       ");
  name.push_back("MOVE WORM         ");
  name.push_back("INSERT KINK       ");
  name.push_back("DELETE KINK       ");
  name.push_back("GLUE WORM         ");
  os << "\n\n# UPDATE STATISTICS";
  os << "\n" << "# col 1 : all updates"
     << "\n" << "# col 2 : impossible updates"
     << "\n" << "# col 3 : rejected updates"
     << "\n" << "# col 4 : accepted updates"
     << "\n" << "# col 5 : acceptance factor with respect to Metropolis ratio only"
     << "\n" << "# col 6 : acceptance factor with respect to all attempts.\n";
  for (size_t i = 0; i < update_count; i++) {
    os << "\n" << name[i]
       << "\t" << update_statistics[total_attempted][i]
       << "\t" << update_statistics[impossible][i]
       << "\t" << update_statistics[rejected][i]
       << "\t" << update_statistics[accepted][i]
       << "\t" << update_statistics[accepted][i] / (update_statistics[total_attempted][i] - update_statistics[impossible][i])
       << "\t" << update_statistics[accepted][i] / update_statistics[total_attempted][i];
  }
  os << "\n\n";
  double s = 0;
  for (size_t i = 0; i < update_count; i++) {
    s += update_statistics[total_attempted][i];
  }
  os << "# Total number of steps : " << s << "\t or 10^" << log10(s) << "\n";
  
}

