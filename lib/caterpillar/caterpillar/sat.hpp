/*------------------------------------------------------------------------------
| This file is distributed under the MIT License.
| See accompanying file /LICENSE for details.
| Author(s): Giulia Meuli
*-----------------------------------------------------------------------------*/
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
  pebble_solver( Network const& net, uint8_t pebbles )
      : _net( net ),
        _pebbles( 4 ),
        _nr_gates( net.num_gates())
  {

    net.foreach_gate( [&]( auto a, auto i ) {
    
      gate_to_index[a] = i;
      index_to_gate[i] = a;
      
    });

    net.foreach_po( [&]( auto po ) {
      o_set.push_back( net.get_node( po ) );
    });

    for(auto e : gate_to_index)
      std::cout <<"gate: " << e.first << ", index: " << gate_to_index[e.first] << std::endl;

    extra = _pebbles*(_nr_gates - _pebbles);
  }

  inline void add_edge_clause( int const& p, int const& p_n, int const& ch, int const& ch_n )
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
    int lit[1];

    solver.set_nr_vars(_nr_gates + extra);
    /* set constraint that everything is unpebbled */
    for ( int v = 0; v < _nr_gates; v++ )
    {
      lit[0] = pabc::Abc_Var2Lit( v, 1 );// zero is not negated
      solver.add_clause( lit, lit + 1 );
    }

    std::cout << "initialize" << std::endl;
    
  }

  void add_step()
  {
    _nr_steps++;
    solver.set_nr_vars( (_nr_gates + extra) * (1 + _nr_steps ) );
    std::cout << "set " << (_nr_gates + extra) * (1 + _nr_steps ) << " vars" << std::endl;

    /* encode move */
    _net.foreach_gate( [&]( auto n, auto i ) //n node i counter
    {
      auto p = ( _nr_steps - 1 ) * (_nr_gates + extra) + i;
      auto p_next = (_nr_steps)*(_nr_gates + extra) + i;

      _net.foreach_fanin( n, [&]( auto ch ) {
        
        auto ch_node = _net.get_node( ch );
        if ( !_net.is_pi( ch_node ) )
        {
          auto ch = gate_to_index[ch_node] + ( _nr_steps - 1 ) * (_nr_gates + extra);
          auto ch_next = gate_to_index[ch_node] + (_nr_steps)* (_nr_gates + extra);
          add_edge_clause( p, p_next, ch, ch_next );
        }
      });
    });

    std::cout << "added step" << std::endl;
    /* cardinality constraint */

    if(_nr_gates > _pebbles)
    {
      /* var declaration */
      std::vector<std::vector<int>> card_vars (_nr_gates - _pebbles);
      auto id_start = _nr_steps*(_nr_gates+extra) + _nr_gates; 
      for(int j = 0 ; j< (_nr_gates - _pebbles) ; j++)
      {
        for(int k = 0 ; k< _pebbles; k++ )
        {
          card_vars[j].push_back(id_start);
          id_start++;
        }
      }

      std::cout << "declared var" << std::endl;

      /* constraint */
      //..
      for (int j = 0; j < _nr_gates - _pebbles - 1; j++)
      {
        for (int k = 0; k < _pebbles; k++)
        {
          int to_or[2];
          to_or[0] = pabc::Abc_Var2Lit( card_vars[j][k], 1);
          to_or[1] = pabc::Abc_Var2Lit( card_vars[j+1][k], 0);
          solver.add_clause(to_or, to_or+2);

          //std::cout << "added clause 1: " << card_vars[j][k] << " || " << card_vars[j+1][k] << std::endl;
        }
      }

      //..
      for (int j = 0; j < _nr_gates - _pebbles; j++)
      {

        for(int kp = 0; kp <= _pebbles ; kp++)
        {
          int k = kp-1;
          std::cout << "k: "<< k << " j: " << j << std::endl;
          int to_var_or[3];

          if(k == -1)
          {
            to_var_or[0] = pabc::Abc_Var2Lit( _nr_steps*(_nr_gates + extra) + j + k + 1, 1);
            to_var_or[1] = pabc::Abc_Var2Lit( card_vars[j][k+1], 0);
            solver.add_clause(to_var_or, to_var_or+2);

            std::cout << "added clause 1: " << _nr_steps*(_nr_gates + extra) + j + k + 1 << " || " << card_vars[j][k+1] << std::endl;
          }
          else if(k == _pebbles - 1)
          {
            to_var_or[0] = pabc::Abc_Var2Lit( _nr_steps*(_nr_gates + extra) + j + k + 1 , 1);
            to_var_or[1] = pabc::Abc_Var2Lit( card_vars[j][k], 1);
            solver.add_clause(to_var_or, to_var_or+2);
            std::cout << "added clause 1: " << _nr_steps*(_nr_gates + extra) + j + k + 1 << " || " << card_vars[j][k] << std::endl;

          }
          else 
          {
            to_var_or[0] = pabc::Abc_Var2Lit( _nr_steps*(_nr_gates + extra) + j + k + 1, 1);
            to_var_or[1] = pabc::Abc_Var2Lit( card_vars[j][k], 1);
            to_var_or[2] = pabc::Abc_Var2Lit( card_vars[j][k+1], 0);
            solver.add_clause(to_var_or, to_var_or+3);

            std::cout << "added clause 1: " << _nr_steps*(_nr_gates + extra) + j + k + 1 << " || " <<card_vars[j][k] << " || " << card_vars[j][k+1] << std::endl;

          }
        }
      }
      
    }

    std::cout << "constraint2 DONE" << std::endl;
  }

  percy::synth_result solve()
  {
    std::vector<int> p (_nr_gates);
    _net.foreach_gate( [&]( auto n, auto i ) //n node i counter
    {
      if ( std::find( o_set.begin(), o_set.end(), n ) != o_set.end() )//is out -> pos state
      {
        p[i] = pabc::Abc_Var2Lit( (_nr_steps)*(_nr_gates + extra) + i, 0 );
      }
      else //is NOT out -> neg state
      {
        p[i] = pabc::Abc_Var2Lit( (_nr_steps)*(_nr_gates + extra) + i, 1 );
      }
    });
    return solver.solve( &p[0], &p[0] + _nr_gates, 0 );
  }
  
  void restart()
  {
    solver.restart();
  }

  void set_nr_vars(int num_steps)
  {
    solver.set_nr_vars(_nr_gates+num_steps*_nr_gates);
  }

  Steps extract_result()
  {
    std::vector<std::vector<int>> vals_step (_nr_steps+1);
    std::vector<int> computed;
    Steps steps;

    
    for ( int i = 0; i <= _nr_steps; i++ )
    {
      for ( int j = 0; j < _nr_gates; j++ )
      {
        vals_step[i].push_back( solver.var_value( (i*(extra+_nr_gates))+j ) );
      }
    }

    std::cout<< std::endl; 
    int s =0;
    for(auto step: vals_step)
    {
      std::cout << "step " << s << ": ";
      s++;
      for (auto node: step)
      {
        std::cout << node;
      }
      std::cout<< std::endl;

    }
  
    for ( int s = 1; s <= _nr_steps; s++ )
    {
      for ( int n = 0; n < _nr_gates; n++ )
      {
        if ( vals_step[s][n] != vals_step[s - 1][n] )
        {
          auto id = std::find( computed.begin(), computed.end(), n );
          bool inplace = false;
          mockturtle::node<Network> target;

          kitty::dynamic_truth_table tt = _net.node_function(gate_to_index[n]);
			    kitty::dynamic_truth_table clone (tt.num_vars());
			    kitty::create_symmetric( clone, mockturtle::detail::odd_bits());
			
          if( (std::find( o_set.begin(), o_set.end(), index_to_gate[n] ) == o_set.end())
              && (clone == tt) )//not out and xor function
          {
            _net.foreach_fanin(index_to_gate[n], [&] (auto ch_signal)
            {
              auto ch_node = _net.get_node(ch_signal);
              if (_net.fanout_size(ch_node) == 1)
              {
                inplace = true;
                target =  ch_node ;
              }
            });
          }
            
          if ( id != computed.end() )
          {
            if(inplace)
              steps.push_back( {index_to_gate[n], uncompute_inplace_action{static_cast<uint32_t>(
                                               target )}} );
            else 
              steps.push_back( {index_to_gate[n], uncompute_action{}} );
            computed.erase( id );
          }
          else
          {
            if(inplace)
              steps.push_back( {index_to_gate[n], compute_inplace_action{static_cast<uint32_t>(
                                               target )}} );
            else 
              steps.push_back( {index_to_gate[n], compute_action{}} );
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
  Network const& _net;
  uint32_t _pebbles;
  uint32_t _nr_gates;
  uint32_t _nr_steps = 0;
  uint32_t extra;
};

template<typename Network>
class pebble_solver_man
{
public:
  using Steps = std::vector<std::pair<mockturtle::node<Network>, mapping_strategy_action>>;

  pebble_solver_man( Network const& ntk, int const pebbles )
      : _ntk( ntk ), _pebbles( pebbles ) {}

  bool solve_nr_steps(pebble_solver<Network>& p_solver)
  { 
    if(p_solver.solve() == percy::success)
    {
      std::cout << "found solution" << std::endl;
      steps = p_solver.extract_result();
      return false;
    }
    else
    {
      return true;
    }
  }

  Steps get_steps()
  {
    pebble_solver<Network> p_solver( _ntk, _pebbles );

    p_solver.initialize(); // initial clauses and set the number of variables

    while(solve_nr_steps(p_solver))
    {
      p_solver.add_step();
    }

    return steps;
  }

private:
  Network const& _ntk;
  int _pebbles;

  Steps steps;
};

} // namespace caterpillar
