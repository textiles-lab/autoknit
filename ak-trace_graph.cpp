#include "pipeline.hpp"

#include <iostream>

void ak::trace_graph(
	ak::RowColGraph const &graph, //in: row-column graph
	std::vector< ak::TracedStitch > *traced_, //out:traced list of stitches
	ak::Model *DEBUG_model_ //in (optional): model
) {
	std::vector< ak::RowColGraph::Vertex > const &vertices = graph.vertices;

	assert(traced_);
	auto &traced = *traced_;
	traced.clear();

	//PARANOIA:
	for (auto const &v : vertices) {
		uint32_t vi = &v - &vertices[0];
		if (v.row_in != -1U) {
			assert(v.row_in < vertices.size());
			assert(vertices[v.row_in].row_out == vi);
		}
		if (v.row_out != -1U) {
			assert(v.row_out < vertices.size());
			assert(vertices[v.row_out].row_in == vi);
		}
		if (v.col_in[0] != -1U) {
			assert(v.col_in[0] < vertices.size());
			assert(v.col_in[0] != v.col_in[1]);
			assert(vertices[v.col_in[0]].col_out[0] == vi || vertices[v.col_in[0]].col_out[1] == vi);
		}
		if (v.col_in[1] != -1U) {
			assert(v.col_in[1] < vertices.size());
			assert(v.col_in[1] != v.col_in[0]);
			assert(vertices[v.col_in[1]].col_out[0] == vi || vertices[v.col_in[1]].col_out[1] == vi);
		}
		if (v.col_out[0] != -1U) {
			assert(v.col_out[0] < vertices.size());
			assert(v.col_out[0] != v.col_out[1]);
			assert(vertices[v.col_out[0]].col_in[0] == vi || vertices[v.col_out[0]].col_in[1] == vi);
		}
		if (v.col_out[1] != -1U) {
			assert(v.col_out[1] < vertices.size());
			assert(v.col_out[1] != v.col_out[0]);
			assert(vertices[v.col_out[1]].col_in[0] == vi || vertices[v.col_out[1]].col_in[1] == vi);
		}
	}
	//end PARANOIA

	//--------------------

	struct VertexInfo {
		uint32_t row = -1U;
		uint32_t knits = 0;
		uint32_t last_stitch = -1U;
	};

	std::vector< VertexInfo > info(vertices.size());
	std::vector< uint32_t > row_pending;

	//divide vertices into rows:
	for (uint32_t seed = 0; seed < vertices.size(); ++seed) {
		if (info[seed].row != -1U) continue;

		//make a new course by making a new slot in the pending array:
		uint32_t row = row_pending.size();
		row_pending.emplace_back(0);

		std::vector< uint32_t > todo;
		info[seed].row = row;
		todo.emplace_back(seed);

		while (!todo.empty()) {
			uint32_t at = todo.back();
			todo.pop_back();
			assert(info[at].row == row);
			for (uint32_t n : {vertices[at].row_in, vertices[at].row_out}) {
				if (n != -1U) {
					if (info[n].row != row) {
						assert(info[n].row == -1U);
						info[n].row = row;
						todo.emplace_back(n);
					}
				}
			}
		}
	}

	std::cout << "Found " << row_pending.size() << " rows." << std::endl;

	//count how many stitches rows are waiting on:
	for (auto const &v : vertices) {
		for (auto n : v.col_out) {
			if (n != -1U) {
				row_pending[info[n].row] += 1;
			}
		}
	}

	//------ actual tracing --------

	constexpr ak::TracedStitch::Dir Forward = ak::TracedStitch::CCW;
	constexpr ak::TracedStitch::Dir Backward = ak::TracedStitch::CW;

	uint32_t fresh_yarn_id = 0;

	auto trace_yarn = [&]() -> bool {
		uint32_t at = -1U;
		ak::TracedStitch::Dir dir = Forward;

		uint32_t prev_stitch = -1U;
		uint32_t yarn = fresh_yarn_id++;
	
		//move to 'next' and make a stitch there:
		auto make_stitch = [&](uint32_t next, ak::TracedStitch::Type type) {

			//DEBUG:
			std::cout << "Make " << char(type) << " at " << next; std::cout.flush(); //DEBUG

			at = next;

			assert(type != ak::TracedStitch::Knit || row_pending[info[at].row] == 0); //tuck on next row sometimes
			assert(type != ak::TracedStitch::Knit || info[at].knits < 2);

			//some type lawyering:
			ak::TracedStitch::Type fancy_type;
			if (type == ak::TracedStitch::Knit) {
				if (info[at].last_stitch == -1U) {
					if (vertices[at].col_in[0] == -1U && vertices[at].col_in[1] == -1U) {
						fancy_type = ak::TracedStitch::Start;
					} else if (vertices[at].col_in[0] != -1U && vertices[at].col_in[1] != -1U) {
						fancy_type = ak::TracedStitch::Decrease;
					} else {
						fancy_type = type;
					}
				} else {
					assert(info[at].knits <= 1); //could be a miss/tuck before as well
					if (vertices[at].col_out[0] == -1U && vertices[at].col_out[1] == -1U) {
						fancy_type = ak::TracedStitch::End;
					} else if (vertices[at].col_out[0] != -1U && vertices[at].col_out[1] != -1U) {
						fancy_type = ak::TracedStitch::Increase;
					} else {
						fancy_type = type;
					}
				}
			} else { //tuck/miss
				assert(type == ak::TracedStitch::Tuck || type == ak::TracedStitch::Miss);
				if (info[at].last_stitch == -1U) {
					//make sure there was a previous stitch that was an increase:
					assert(vertices[at].col_in[0] != -1U && vertices[at].col_in[1] == -1U);
					assert(vertices[at].col_in[0] < vertices.size());
					assert(vertices[vertices[at].col_in[0]].col_out[0] != -1U && vertices[vertices[at].col_in[0]].col_out[1] != -1U);
					assert(vertices[at].col_in[0] < info.size());
					assert(info[vertices[at].col_in[0]].last_stitch != -1U);
				} else {
					ak::TracedStitch::Type last_type = traced[info[at].last_stitch].type;
					assert(last_type != ak::TracedStitch::Increase);
				}
				fancy_type = type;
			}
			std::cout << " (became " << char(fancy_type) << ")" << std::endl; //DEBUG

			//build stitch:
			TracedStitch ts;
			ts.yarn = yarn;
			if (info[at].last_stitch != -1U) {
				ts.ins[0] = info[at].last_stitch;
			} else {
				if (vertices[at].col_in[0] != -1U) ts.ins[0] = info[vertices[at].col_in[0]].last_stitch;
				if (vertices[at].col_in[1] != -1U) ts.ins[1] = info[vertices[at].col_in[1]].last_stitch;
				if (dir == Backward && ts.ins[0] != -1U && ts.ins[1] != -1U) {
					std::swap(ts.ins[0], ts.ins[1]);
				}
			}
			ts.type = fancy_type;
			ts.dir = dir;
			ts.vertex = at;


			//update info:
			if (type == ak::TracedStitch::Knit) {
				info[at].knits += 1;
				if (info[at].knits == 2) {
					for (auto n : vertices[at].col_out) {
						if (n != -1U) {
							assert(row_pending[info[n].row] > 0);
							row_pending[info[n].row] -= 1;
						}
					}
				}
			}
			info[at].last_stitch = traced.size();

			//store stitch:
			prev_stitch = traced.size();
			traced.emplace_back(ts);
		};

		auto knit = [&](uint32_t next) { make_stitch(next, ak::TracedStitch::Knit); };
		auto tuck = [&](uint32_t next) { make_stitch(next, ak::TracedStitch::Tuck); };
		auto miss = [&](uint32_t next) { make_stitch(next, ak::TracedStitch::Miss); };

		//Rule 1: start by knitting a ready but not knit-twice node:
		{
			uint32_t found = -1U;
			for (uint32_t vi = 0; vi < info.size(); ++vi) {
				if (row_pending[info[vi].row] == 0 && info[vi].knits < 2) {
					found = vi;
					break;
				}
			}
			if (found == -1U) return false;
			//now shove 'found' to one end of a chain of ready stitches:
			auto adv = [&](uint32_t *vi) {
				uint32_t prev = vertices[*vi].row_in;
				//NOTE: probably some special cases to consider here around the ends of short rows; will ignore them for now.

				//if previous stitch exists but is already knit on, consider moving to a child of that stitch:
				while (prev != -1U && info[prev].knits >= 2) {
					assert(row_pending[info[prev].row] == 0); //it's the same row, so... yeah.
					uint32_t found = -1U;
					for (uint32_t o : {1 , 0}) {
						uint32_t child = vertices[prev].col_out[o];
						if (child == -1U) continue;
						if (row_pending[info[child].row] != 0) continue;
						found = child;
						break;
					}
					prev = found;
				}
				//no previous stitch (or previous stitch child) was found that could support knits:
				if (prev == -1U) return false;
				//previous stitch found that could be knit:
				assert(row_pending[info[prev].row] == 0);
				assert(info[prev].knits < 2);
				*vi = prev;
				return true;
			};

			uint32_t found2 = found;
			while (true) {
				if (!adv(&found) || found == found2) break;
				if (!adv(&found) || found == found2) break;
				bool ret = adv(&found2);
				assert(ret);
				if (found == found2) break;
			}

			//Swap direction based on row-wise neighbors:
			if (vertices[found].row_out == -1U || info[vertices[found].row_out].knits == 2) dir = Backward;
			else dir = Forward;
			knit(found);
		}
		if (at == -1U) return false;

		auto get_next = [&](uint32_t v) {
			assert(v < vertices.size());
			return (dir == Forward ? vertices[v].row_out : vertices[v].row_in);
		};

		/* unused
		auto get_prev = [&](uint32_t v) {
			assert(v < vertices.size());
			return (dir == Forward ? vertices[v].row_in : vertices[v].row_out);
		};
		*/

		auto get_prev_child = [&](uint32_t v) {
			assert(v < vertices.size());
			return (dir == Forward ? vertices[v].col_out[0] : vertices[v].col_out[1]);
		};

		auto get_next_child = [&](uint32_t v) {
			assert(v < vertices.size());
			return (dir == Forward ? vertices[v].col_out[1] : vertices[v].col_out[0]);
		};

		auto get_prev_parent = [&](uint32_t v) {
			assert(v < vertices.size());
			return (dir == Forward ? vertices[v].col_in[0] : vertices[v].col_in[1]);
		};

		auto get_next_parent = [&](uint32_t v) {
			assert(v < vertices.size());
			return (dir == Forward ? vertices[v].col_in[1] : vertices[v].col_in[0]);
		};


		auto is_covered = [&](uint32_t v) {
			//'covered' == vertex v already has stitches on successors
			assert(v < vertices.size());
			for (auto n : vertices[v].col_out) {
				if (n != -1U && info[n].knits != 0) return true;
			}
			return false;
		};


		//Rule 2: move to next row if next row is ready
		auto rule2 = [&]() -> bool {
			uint32_t up = -1U;
			for (uint32_t n : {get_next_child(at), get_prev_child(at)}) {
				if (n != -1U && row_pending[info[n].row] == 0) {
					up = n;
					break;
				}
			}
			if (up == -1U) return false;
			//'at' is below 'up' which is ready

			{ //if 'up' has a next neighbor, just knit over to it:
				uint32_t up_next = get_next(up);
				if (up_next != -1U) {
					knit(up_next);
					return true;
				}
			}

			//otherwise, need to (maybe) tuck, then turn:
			uint32_t next = get_next(at);

			if (next != -1U && is_covered(next)) {
				next = -1U;
				std::cout << "NOTE: not tucking because neighbor is covered." << std::endl;
			}

			if (next != -1U && info[next].knits == 2 && vertices[next].col_out[0] != -1U && vertices[next].col_out[1] != -1U) {
				next = get_prev_child(next);
				std::cout << "NOTE: tucking on child of next because of increase." << std::endl;
			}

			if (next != -1U && info[next].last_stitch != -1U && traced[info[next].last_stitch].type == ak::TracedStitch::End) {
				std::cout << "NOTE: not tucking on next because it is an end." << std::endl;
				next = -1U;
			}


			if (next != -1U) {
				std::cout << "  TUCKING[2] at " << next << " which has " << info[next].knits << " knits." << std::endl;
			}

			//tuck 'next', turn, knit 'up':
			if (next != -1U) tuck(next);
			dir = (dir == Forward ? Backward : Forward);
			if (next != -1U) miss(next);
			knit(up);
			return true;
		};

		//Rule 3: continue along row if no column edge to ready next stitch:
		auto rule3 = [&]() -> bool {
			//is there a next stitch in this row that needs knitting?
			uint32_t next = get_next(at);
			if (next == -1U || info[next].knits >= 2) return false;

			//yep, so knit it:
			knit(next);
			return true;
		};

		//Rule 4: tuck and turn at the end of short rows:
		auto rule4 = [&]() -> bool {
			//must be at a stitch knit only once:
			if (info[at].knits == 2) return false;
			assert(info[at].knits == 1);

			//must be no next stitch in this row:
			uint32_t next = get_next(at);
			if (next != -1U) return false;

			//find parent stitch:
			uint32_t down = -1U;
			for (uint32_t n : {get_next_parent(at), get_prev_parent(at)}) {
				if (n != -1U) {
					down = n;
					break;
				}
			}
			if (down == -1U) return false;
			assert(info[down].knits == 2); //can't be here if parent wasn't knit twice

			//if down has a next neighbor, tuck before turning:
			uint32_t down_next = get_next(down);

			if (down_next != -1U && is_covered(down_next)) {
				std::cout << "NOTE: not tucking in rule4 because down_next is covered." << std::endl;
				down_next = -1U;
			}

			if (down_next != -1U && info[down_next].knits == 2 && vertices[down_next].col_out[0] != -1U && vertices[down_next].col_out[1] != -1U) {
				down_next = get_prev_child(down_next);
				std::cout << "NOTE: tucking on child of down_next because of increase." << std::endl;
			}
			if (down_next != -1U && info[down_next].last_stitch != -1U && traced[info[down_next].last_stitch].type == ak::TracedStitch::End) {
				std::cout << "NOTE: not tucking on down_next because it is an end." << std::endl;
				down_next = -1U;
			}

			if (down_next != -1U) {
				std::cout << "  TUCKING[4] at " << down_next << " which has " << info[down_next].knits << " knits and outs " << int32_t(vertices[down_next].col_out[0]) << " and " << int32_t(vertices[down_next].col_out[1]) << std::endl;
			}


		
			uint32_t here = at; //because tuck() / miss() will change 'at'
			if (down_next != -1U) tuck(down_next);
			dir = (dir == Forward ? Backward : Forward);
			if (down_next != -1U) miss(down_next);
			knit(here);
			return true;
		};

		//Rule 5: walk off the end of short rows:
		auto rule5 = [&]() -> bool {
			//must be at a stitch knit twice:
			if (info[at].knits != 2) return false;

			//find next-most parent stitch with a next neighbor:
			uint32_t par = at;
			uint32_t par_next = -1U;
			while (par_next == -1U) {
				bool found = false;
				for (uint32_t n : {get_next_parent(par), get_prev_parent(par)}) {
					if (n != -1U) {
						par = n;
						found = true;
						break;
					}
				}
				if (!found) return false; //ran out of parents
				par_next = get_next(par);
			}
			assert(par_next != -1U);

			if (info[par_next].knits == 2) return false;
			knit(par_next);

			return true;
		};

		//continue running rules until none fire:
		while (
			   rule2()
			|| rule3()
			|| rule4()
			|| rule5()
		) { /* spin */ }

		//rule6: when all the other rules stop working, end the yarn.

		return true;
	};


	//trace yarns until no more yarns:
	while (trace_yarn()) { /*spin */ }

	//fix up the 'out' pointers from the 'in' pointers:
	auto add_out = [&traced](uint32_t ti, uint32_t b) {
		assert(ti < traced.size());
		if (traced[ti].outs[0] == -1U) {
			traced[ti].outs[0] = b;
		} else if (traced[ti].outs[1] == -1U) {
			traced[ti].outs[1] = b;
		} else {
			assert(traced[ti].outs[0] == -1U || traced[ti].outs[1] == -1U); //gotta have some room!
		}
	};
	for (auto const &ts : traced) {
		uint32_t ti = &ts - &traced[0];
		if (ts.ins[0] != -1U) add_out(ts.ins[0], ti);
		if (ts.ins[1] != -1U) add_out(ts.ins[1], ti);
	}
	
	//sort outs using the handy 'vertex' field:
	for (auto &ts : traced) {
		if (ts.outs[0] != -1U && ts.outs[1] != -1U) {
			if (vertices[ts.vertex].col_out[0] == traced[ts.outs[1]].vertex
			 && vertices[ts.vertex].col_out[1] == traced[ts.outs[0]].vertex) {
				std::swap(ts.outs[0], ts.outs[1]);
			}
			assert( vertices[ts.vertex].col_out[0] == traced[ts.outs[0]].vertex
			     && vertices[ts.vertex].col_out[1] == traced[ts.outs[1]].vertex );

			if (ts.dir == Backward) {
				std::swap(ts.outs[0], ts.outs[1]);
			}
		}
	}


	//set stitch 'at' using model vertex positions:
	if (DEBUG_model_) {
		std::vector< glm::vec3 > at; at.reserve(vertices.size());
		for (auto const &v : vertices) {
			at.emplace_back(v.at.interpolate(DEBUG_model_->vertices));
		}
		assert(at.size() == vertices.size());
		std::vector< glm::vec3 > up; up.reserve(vertices.size());
		for (auto const &v : vertices) {
			glm::vec3 const &v_at = at[&v - &vertices[0]];

			glm::vec3 acc = glm::vec3(0.0f);
			uint32_t sum = 0;
			for (uint32_t i = 0; i < 2; ++i) {
				if (v.col_in[i] != -1U) {
					acc += glm::normalize(v_at - at[v.col_in[i]]);
					sum += 1;
				}
			}
			for (uint32_t i = 0; i < 2; ++i) {
				if (v.col_out[i] != -1U) {
					acc += glm::normalize(at[v.col_out[i]] - v_at);
					sum += 1;
				}
			}
			if (sum == 0) acc = glm::vec3(0.0f, 0.0f, 1.0f);
			else acc = acc / float(sum);

			up.emplace_back(acc);
		}
		assert(up.size() == vertices.size());

		std::vector< uint32_t > count(vertices.size(), 0);
		std::vector< uint32_t > index; index.reserve(traced.size());
		for (auto const &ts : traced) {
			index.emplace_back(count[ts.vertex]);
			count[ts.vertex] += 1;
		}
		assert(index.size() == traced.size());

		for (auto &ts : traced) {
			uint32_t i = index[&ts - &traced[0]];
			float amt = float(i + 0.5f) / float(count[ts.vertex]);
			amt = 0.1f * (2.0f * (amt - 0.5f));
			ts.at = at[ts.vertex] + amt * up[ts.vertex];
		}
	}
}
