/*------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Giulia Meuli
*-----------------------------------------------------------------------------*/
#include <cstdint>
#include <unordered_map>
#include <variant>
#include <vector>

#include <mockturtle/traits.hpp>
#include <mockturtle/utils/node_map.hpp>
#include <percy/solvers/bsat2.hpp>

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
  pebble_solver( Network const& net, uint32_t pebbles )
      : index_to_gate( net.num_gates() ),
        gate_to_index( net ),
        _net( net ),
        _pebbles( pebbles ),
        _nr_gates( net.num_gates() )
  {
    net.foreach_gate( [&]( auto a, auto i ) {
      gate_to_index[a] = i;
      index_to_gate[i] = a;
    } );

    net.foreach_po( [&]( auto po ) {
      o_set.insert( net.get_node( po ) );
    } );

    extra = ( _pebbles < _nr_gates ) ? _pebbles * ( _nr_gates - _pebbles ) : 0;
  }

  inline void add_edge_clause( int p, int p_n, int ch, int ch_n )
  {
    int h[3];
    h[0] = pabc::Abc_Var2Lit( p, 1 );
    h[1] = pabc::Abc_Var2Lit( p_n, 0 );
    h[2] = pabc::Abc_Var2Lit( ch, 0 );
    solver.add_clause( h, h + 3 );

    h[0] = pabc::Abc_Var2Lit( p, 1 );
    h[1] = pabc::Abc_Var2Lit( p_n, 0 );
    h[2] = pabc::Abc_Var2Lit( ch_n, 0 );
    solver.add_clause( h, h + 3 );

    h[0] = pabc::Abc_Var2Lit( p, 0 );
    h[1] = pabc::Abc_Var2Lit( p_n, 1 );
    h[2] = pabc::Abc_Var2Lit( ch, 0 );
    solver.add_clause( h, h + 3 );

    h[0] = pabc::Abc_Var2Lit( p, 0 );
    h[1] = pabc::Abc_Var2Lit( p_n, 1 );
    h[2] = pabc::Abc_Var2Lit( ch_n, 0 );
    solver.add_clause( h, h + 3 );
  }

  void initialize()
  {
    solver.set_nr_vars( _nr_gates + extra );

    /* set constraint that everything is unpebbled */
    for ( int v = 0; v < _nr_gates; v++ )
    {
      int lit = pabc::Abc_Var2Lit( v, 1 ); // zero is not negated
      solver.add_clause( &lit, &lit + 1 );
    }
  }

  void add_step()
  {
    _nr_steps++;
    solver.set_nr_vars( ( _nr_gates + extra ) * ( 1 + _nr_steps ) );

    /* encode move */
    _net.foreach_gate( [&]( auto n, auto i ) {
      auto p = pebble_var( _nr_steps - 1, i );
      auto p_next = pebble_var( _nr_steps, i );

      _net.foreach_fanin( n, [&]( auto ch ) {
        auto ch_node = _net.get_node( ch );
        if ( !_net.is_constant( ch_node ) && !_net.is_pi( ch_node ) )
        {
          const auto ch = pebble_var( _nr_steps - 1, gate_to_index[ch_node] );
          const auto ch_next = pebble_var( _nr_steps, gate_to_index[ch_node] );
          add_edge_clause( p, p_next, ch, ch_next );
        }
      } );
    } );

    /* cardinality constraint */

    if ( ( _pebbles > 0 ) && ( _nr_gates > _pebbles ) )
    {
      /* var declaration */
      std::vector<std::vector<int>> card_vars( _nr_gates - _pebbles );
      auto id_start = _nr_steps * ( _nr_gates + extra ) + _nr_gates;
      for ( int j = 0; j < ( _nr_gates - _pebbles ); j++ )
      {
        for ( int k = 0; k < _pebbles; k++ )
        {
          card_vars[j].push_back( id_start );
          id_start++;
        }
      }

      /* constraint */
      //..
      for ( int j = 0; j < _nr_gates - _pebbles - 1; j++ )
      {
        for ( int k = 0; k < _pebbles; k++ )
        {
          int to_or[2];
          to_or[0] = pabc::Abc_Var2Lit( card_vars[j][k], 1 );
          to_or[1] = pabc::Abc_Var2Lit( card_vars[j + 1][k], 0 );
          solver.add_clause( to_or, to_or + 2 );
        }
      }

      //..
      for ( int j = 0; j < _nr_gates - _pebbles; j++ )
      {

        for ( int kp = 0; kp <= _pebbles; kp++ )
        {
          int k = kp - 1;
          int to_var_or[3];

          if ( k == -1 )
          {
            to_var_or[0] = pabc::Abc_Var2Lit( _nr_steps * ( _nr_gates + extra ) + j + k + 1, 1 );
            to_var_or[1] = pabc::Abc_Var2Lit( card_vars[j][k + 1], 0 );
            solver.add_clause( to_var_or, to_var_or + 2 );
          }
          else if ( k == _pebbles - 1 )
          {
            to_var_or[0] = pabc::Abc_Var2Lit( _nr_steps * ( _nr_gates + extra ) + j + k + 1, 1 );
            to_var_or[1] = pabc::Abc_Var2Lit( card_vars[j][k], 1 );
            solver.add_clause( to_var_or, to_var_or + 2 );
          }
          else
          {
            to_var_or[0] = pabc::Abc_Var2Lit( _nr_steps * ( _nr_gates + extra ) + j + k + 1, 1 );
            to_var_or[1] = pabc::Abc_Var2Lit( card_vars[j][k], 1 );
            to_var_or[2] = pabc::Abc_Var2Lit( card_vars[j][k + 1], 0 );
            solver.add_clause( to_var_or, to_var_or + 3 );
          }
        }
      }
    }
  }

  percy::synth_result solve()
  {
    std::vector<int> p( _nr_gates );
    _net.foreach_gate( [&]( auto n, auto i ) {
      p[i] = pabc::Abc_Var2Lit( pebble_var( _nr_steps, i ), o_set.count( n ) ? 0 : 1 );
    } );
    return solver.solve( &p[0], &p[0] + _nr_gates, 0 );
  }

  inline int pebble_var( int step, int gate )
  {
    return step * ( extra + _nr_gates ) + gate;
  }

  Steps extract_result()
  {
    std::vector<std::vector<int>> vals_step( _nr_steps + 1 );
    Steps steps;

    for ( int i = 0; i <= _nr_steps; i++ )
    {
      for ( int j = 0; j < _nr_gates; j++ )
      {
        const auto value = solver.var_value( pebble_var( i, j ) );
        vals_step[i].push_back( value );
      }
    }

    for ( int s = 1; s <= _nr_steps; s++ )
    {
      auto it = steps.end();

      for ( int n = 0; n < _nr_gates; n++ )
      {
        if ( vals_step[s][n] != vals_step[s - 1][n] )
        {
          bool inplace = false;
          uint32_t target;

          kitty::dynamic_truth_table tt = _net.node_function( gate_to_index[n] );
          auto clone = tt.construct();
          kitty::create_parity( clone );

#if 0
          if ( !o_set.count( index_to_gate[n] ) && clone == tt )
          {
            _net.foreach_fanin(index_to_gate[n], [&] (auto ch_signal)
            {
              auto ch_node = _net.get_node(ch_signal);
              if (_net.fanout_size(ch_node) == 1)
              {
                inplace = true;
                target = _net.node_to_index( ch_node );
                return false;
              }
              return true;
            });
          }
#endif

          if ( !vals_step[s][n] )
          {
            if ( inplace )
            {
              it = steps.insert( it, {index_to_gate[n], uncompute_inplace_action{target}} );
              ++it;
            }
            else
            {
              it = steps.insert( it, {index_to_gate[n], uncompute_action{}} );
              ++it;
            }
          }
          else
          {
            if ( inplace )
            {
              it = steps.insert( it, {index_to_gate[n], compute_inplace_action{target}} );
            }
            else
            {
              it = steps.insert( it, {index_to_gate[n], compute_action{}} );
            }
          }
        }
      }
    }

    return steps;
  }

private:
  std::vector<mockturtle::node<Network>> index_to_gate;
  mockturtle::node_map<int, Network> gate_to_index;
  std::unordered_set<mockturtle::node<Network>> o_set;

  percy::bsat_wrapper solver;
  Network const& _net;
  uint32_t _pebbles;
  uint32_t _nr_gates;
  uint32_t _nr_steps = 0;
  uint32_t extra;
};

} // namespace caterpillar
