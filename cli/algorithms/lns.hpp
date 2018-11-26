#include <alice/alice.hpp>

#include <caterpillar/lhrs.hpp>
#include <tweedledum/algorithms/synthesis/single_target_gates.hpp>

#include "../utils/cirkit_command.hpp"

namespace alice
{

class lns_command : public cirkit::cirkit_command<lns_command, aig_t, mig_t, xag_t, xmg_t, klut_t>
{
public:
  lns_command( environment::ptr& env ) : cirkit::cirkit_command<lns_command, aig_t, mig_t, xag_t, xmg_t, klut_t>( env, "Logic network based hierarchical synthesis", "hierarchical synthesis from {0}" )
  {
    add_option( "--qmapping", qmapping, "qubit mapping", true )->set_type_name( "strategy in {bennett=0, bennett inplace=1, pebbling=2}" );
    add_option( "--gmapping", gmapping, "STG gate mapping", true )->set_type_name( "strategy in {PPRM=0, PKRM=1, spectrum=2}" );
    add_option( "--pebble_limit", pebble_limit, "Maximum number of pebbles for strategy 2" );
    add_option( "--conflict_limit", conflict_limit, "Conflict limit for SAT solver in strategy 2" );
    add_flag("--increment_on_timeout", "increment on timeout");
    add_flag("--decrement_on_success", "decrement on success");
    add_flag( "-p,--progress", "show progress" );
    add_flag( "-v,--verbose", "be verbose" );
    
  }

  template<class Store>
  inline void execute_store()
  {
    /* dispatch by qmapping */
    using LogicNetwork = typename Store::element_type;
    switch ( qmapping )
    {
    default:
      env->err() << "[e] invalid qmapping\n";
      break;
    case 0u:
      execute_store_with_mapping<Store, typename caterpillar::bennett_mapping_strategy<LogicNetwork>>();
      break;
    case 1u:
      execute_store_with_mapping<Store, typename caterpillar::bennett_inplace_mapping_strategy<LogicNetwork>>();
      break;
    case 2u:
      execute_store_with_mapping<Store, typename caterpillar::pebbling_mapping_strategy<LogicNetwork>>();
      break;
    }
  }

private:
  template<class Store, class MappingStrategy>
  inline void execute_store_with_mapping()
  {
    using LogicNetwork = typename Store::element_type;

    ps.verbose = is_set( "verbose" );
    if ( is_set( "pebble_limit" ) )
    {
      ps.mapping_ps.pebble_limit = pebble_limit;
    }
    if ( is_set( "conflict_limit" ) )
    {
      ps.mapping_ps.conflict_limit = conflict_limit;
    }

    ps.mapping_ps.increment_on_timeout = is_set("increment_on_timeout");
    ps.mapping_ps.decrement_on_success = is_set("decrement_on_success");

    ps.mapping_ps.progress = is_set( "progress" );
    auto& circs = store<qcircuit_t>();
    if ( circs.empty() || is_set( "new" ) )
    {
      circs.extend();
    }
    circs.current() = qcircuit_t();

    /* dispatch by gmapping */
    switch ( gmapping )
    {
    default:
      env->err() << "[e] invalid gmapping\n";
      break;
    case 0u:
      success = caterpillar::logic_network_synthesis<qcircuit_t, LogicNetwork, MappingStrategy>( circs.current(), *( store<Store>().current() ), tweedledum::stg_from_pprm(), ps, &st );
      break;
    case 1u:
      success = caterpillar::logic_network_synthesis<qcircuit_t, LogicNetwork, MappingStrategy>( circs.current(), *( store<Store>().current() ), tweedledum::stg_from_pkrm(), ps, &st );
      break;
    case 2u: {
      tweedledum::stg_from_spectrum_params stg_ps;
      stg_ps.lin_comb_synth_behavior = tweedledum::stg_from_spectrum_params::never;
      success = caterpillar::logic_network_synthesis<qcircuit_t, LogicNetwork, MappingStrategy>( circs.current(), *( store<Store>().current() ), tweedledum::stg_from_spectrum( stg_ps ), ps, &st );
    } break;
    }
  }

public:
  nlohmann::json log() const override
  {
    return {
        {"qmapping", qmapping},
        {"gmapping", gmapping},
        {"success", success},
        {"required_ancillae", st.required_ancillae},
        {"time_total", mockturtle::to_seconds( st.time_total )}};
  }

private:
  caterpillar::logic_network_synthesis_params ps;
  caterpillar::logic_network_synthesis_stats st;

  unsigned qmapping{0u};
  unsigned gmapping{0u};
  uint32_t pebble_limit;
  uint32_t conflict_limit;
  bool success{true};
};

ALICE_ADD_COMMAND( lns, "Synthesis" )

} // namespace alice
