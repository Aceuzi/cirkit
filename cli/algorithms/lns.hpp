#include <alice/alice.hpp>

#include <tweedledum/algorithms/synthesis/logic_network_synthesis.hpp>

#include "../utils/cirkit_command.hpp"

namespace alice
{

class lns_command : public cirkit::cirkit_command<lns_command, mig_t, xmg_t, klut_t>
{
public:
  lns_command( environment::ptr& env ) : cirkit::cirkit_command<lns_command, mig_t, xmg_t, klut_t>( env, "Logic network based hierarchical synthesis", "hierarchical synthesis from {0}" )
  {
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
    logic_network_synthesis( circs.current(), *( store<Store>().current() ), ps );
  }

private:
  tweedledum::logic_network_synthesis_params ps;
};

ALICE_ADD_COMMAND( lns, "Synthesis" )

} // namespace alice
