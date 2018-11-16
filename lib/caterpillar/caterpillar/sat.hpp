#include <percy/solvers/bsat2.hpp>
#include <variant>

namespace caterpillar
{

struct compute_action
{
};
struct uncompute_action
{
};
struct compute_inplace_action
{
  uint32_t target_index;
};
struct uncompute_inplace_action
{
  uint32_t target_index;
};

using mapping_strategy_action = std::variant<compute_action, uncompute_action, compute_inplace_action, uncompute_inplace_action>;

template<typename Network>
class pebble_solver
{
  using Steps = std::vector<std::pair<mockturtle::node<Network>, mapping_strategy_action>>;

public:
  pebble_solver( Network net, uint8_t pebbles )
      : _net( net ),
        _pebbles( pebbles ),
        _nr_steps( 0 ),
        _nr_gates( net.num_gates() )
  {
    _net.foreach_gate( [&]( auto a, auto i ) {
      gate_to_index[a] = i;
      index_to_gate[i] = a;
    } );

    _net.foreach_po( [&]( auto po ) {
      o_set.push_back( _net.get_node( po ) );
    } );
  }

  inline void add_edge_clause( int const& p, int const& p_n, int const& ch, int const& ch_n )
  {
    int h[3];
    h[0] = pabc::Abc_Var2Lit( p, 0 );
    h[1] = pabc::Abc_Var2Lit( p_n, 1 );
    h[2] = pabc::Abc_Var2Lit( ch, 1 );
    solver.add_clause( h, h + 3 );

    h[0] = pabc::Abc_Var2Lit( p, 1 );
    h[1] = pabc::Abc_Var2Lit( p_n, 0 );
    h[2] = pabc::Abc_Var2Lit( ch, 1 );
    solver.add_clause( h, h + 3 );

    h[0] = pabc::Abc_Var2Lit( p, 0 );
    h[1] = pabc::Abc_Var2Lit( p_n, 1 );
    h[2] = pabc::Abc_Var2Lit( ch_n, 1 );
    solver.add_clause( h, h + 3 );

    h[0] = pabc::Abc_Var2Lit( p, 1 );
    h[1] = pabc::Abc_Var2Lit( p_n, 0 );
    h[2] = pabc::Abc_Var2Lit( ch_n, 1 );
    solver.add_clause( h, h + 3 );
  }

  void initialize()
  {
    int lit[1];

    /* set constraint that everything is unpebbled */
    for ( int v = 0; v < _nr_gates; v++ )
    {
      lit[0] = pabc::Abc_Var2Lit( v, 0 );
      solver.add_clause( lit, lit + 1 );
    }
  }

  void add_step()
  {
    _nr_steps++;

    /* encode move */
    _net.foreach_gate( [&]( auto n, auto i ) //n node i counter
                       {
                         auto p = ( _nr_steps - 1 ) * _nr_gates + i;
                         auto p_next = (_nr_steps)*_nr_gates + i;

                         _net.foreach_fanin( n, [&]( auto ch ) {
                           // if pi not use it
                           auto ch_node = _net.get_node( ch );
                           if ( !_net.is_pi( ch_node ) )
                           {
                             auto ch = gate_to_index[ch_node] + ( _nr_steps - 1 ) * _nr_gates;
                             auto ch_next = gate_to_index[ch_node] + (_nr_steps)*_nr_gates;

                             add_edge_clause( p, p_next, ch, ch_next );
                           }
                         } );
                       } );

    /* cardinality constraint */
  }

  void add_final()
  {
    int p[1];
    _net.foreach_gate( [&]( auto n, auto i ) //n node i counter
                       {
                         if ( std::find( o_set.begin(), o_set.end(), n ) != o_set.end() )
                         {
                           p[0] = pabc::Abc_Var2Lit( (_nr_steps)*_net.num_gates() + i, 1 );
                         }
                         else
                         {
                           p[0] = pabc::Abc_Var2Lit( (_nr_steps)*_net.num_gates() + i, 0 );
                         }

                         solver.add_clause( p, p + 1 );
                       } );
  }

  void remove_final()
  {
    int p[_nr_gates];
    _net.foreach_gate( [&]( auto n, auto i ) {
      if ( std::find( o_set.begin(), o_set.end(), n ) != o_set.end() )
      {
        p[i] = pabc::Abc_Var2Lit( (_nr_steps)*_net.num_gates() + i, 0 );
      }
      else
      {
        p[i] = pabc::Abc_Var2Lit( (_nr_steps)*_net.num_gates() + i, 1 );
      }
    } );
    solver.add_clause( p, p + _nr_gates );

  }

  percy::synth_result solve()
  {
    return solver.solve( 0 );
  }

  Steps extract_result()
  {
    std::vector<std::vector<int>> vals_step;
        std::vector<int>
            computed;
    Steps steps;
    int c = 0;
    for ( int i = 0; i < _nr_steps; i++ )
    {
      for ( int j = 0; j < _nr_gates; j++ )
      {
        vals_step[i].push_back( solver.var_value( c ) );
        c++;
      }
    }

    for ( int s = 1; s < _nr_steps; s++ )
    {
      for ( int n = 0; n < _nr_gates; n++ )
      {
        if ( vals_step[s][n] != vals_step[s - 1][n] )
        {
          auto id = std::find( computed.begin(), computed.end(), n );
          if ( id != computed.end() )
          {
            steps.push_back( { n, uncompute_action{} } );
            computed.erase( id );
          }
          else
          {
            steps.push_back( { n, compute_action{} } );
            computed.push_back( n );
          }
        }
      }
    }

    return steps;
  }

private:
  std::map<int, int> index_to_gate;
  std::map<int, int> gate_to_index;
  std::vector<int> o_set;

  percy::bsat_wrapper solver;
  uint32_t _pebbles;
  Network _net;
  uint32_t _nr_steps;
  uint32_t _nr_gates;
};

template<typename Network>
class pebble_solver_man
{
public:
  using Steps = std::vector<std::pair<mockturtle::node<Network>, mapping_strategy_action>>;

  pebble_solver_man( Network const& ntk, int const pebbles )
      : _ntk( ntk ), _pebbles( pebbles ) {}

  Steps get_steps()
  {
    pebble_solver<Network> p_solver( _ntk, _pebbles );

    p_solver.initialize();
    p_solver.add_final();

    while ( p_solver.solve() == percy::failure )
    {
      p_solver.remove_final();
      p_solver.add_step();
      p_solver.add_final();
    }
    return p_solver.extract_result();
  }

private:
  Network _ntk;
  int _pebbles;
};

} // namespace caterpillar
