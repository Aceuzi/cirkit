/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @file write_graphviz.hpp
 *
 * @brief Use graphviz to render and print the graph
 *
 * @author Heinz Riener
 * @since  2.0
 */

#ifndef WRITE_GRAPHVIZ_HPP
#define WRITE_GRAPHVIZ_HPP

#if ADDON_GRAPHVIZ

#include <classical/aig.hpp>
#include <classical/graphviz.hpp>
#include <string>
#include <map>

namespace cirkit
{

void compute_graphviz_layout( gv_graph& gv, const aig_graph& aig, const std::string& layout_algorithm, const std::string& render_format,
                              std::map< aig_node, gv_node >& node_map, std::map< aig_edge, gv_edge >& edge_map );

void write_graphviz( const aig_graph& aig, const std::string& layout_algorithm, const std::string& graphviz_renderer, std::ostream& os );
void write_graphviz( const aig_graph& aig, const std::string& layout_algorithm, const std::string& render_format,
                     std::map< aig_node, gv_node >& node_map, std::map< aig_edge, gv_edge >& edge_map, std::ostream& os );

}

#endif

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
