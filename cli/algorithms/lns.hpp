#include <alice/alice.hpp>

#include <caterpillar/lhrs.hpp>

#include "../utils/cirkit_command.hpp"

namespace alice
{

class lns_command : public cirkit::cirkit_command<lns_command, aig_t, mig_t, xag_t, xmg_t, klut_t>
{
public:
  lns_command( environment::ptr& env ) : cirkit::cirkit_command<lns_command, aig_t, mig_t, xag_t, xmg_t, klut_t>( env, "Logic network based hierarchical synthesis", "hierarchical synthesis from {0}" )
  {
    add_flag( "--outofplace", "use always out-of-place mapping" );
    add_flag( "--pebbling", "use always out-of-place pebbling" );
    add_flag( "-v,--verbose", "be verbose" );
  }

  template<class Store>
  inline void execute_store()
  {
    ps.verbose = is_set( "verbose" );
    auto& circs = store<qcircuit_t>();
    if ( circs.empty() || is_set( "new" ) )
    {
      circs.extend();
    }
    circs.current() = qcircuit_t();

    using LogicNetwork = typename Store::element_type;
    if ( is_set( "pebbling" ) )
    {
      caterpillar::logic_network_synthesis<qcircuit_t, LogicNetwork, typename caterpillar::pebbling_mapping_strategy<LogicNetwork>>( circs.current(), *( store<Store>().current() ), ps );
    }
    if ( is_set( "outofplace" ) )
    {
      caterpillar::logic_network_synthesis<qcircuit_t, LogicNetwork, typename caterpillar::bennett_mapping_strategy<LogicNetwork>>( circs.current(), *( store<Store>().current() ), ps );
    }
    else
    {
      caterpillar::logic_network_synthesis<qcircuit_t, LogicNetwork, typename caterpillar::bennett_inplace_mapping_strategy<LogicNetwork>>( circs.current(), *( store<Store>().current() ), ps );
    }
  }

private:
  caterpillar::logic_network_synthesis_params ps;
};

ALICE_ADD_COMMAND( lns, "Synthesis" )

} // namespace alice
