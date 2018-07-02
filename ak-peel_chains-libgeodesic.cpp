#include "pipeline.hpp"

extern "C" {
#include <hmContext.h>
#include <hmTriMesh.h>
#include <hmTriDistance.h>
#include <hmUtility.h>
}

#include <iostream>

void ak::peel_chains(
	ak::Parameters const &parameters,
	ak::Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< ak::EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< ak::EmbeddedVertex > > *next_chains_ //out: next chains (may be different size than active_chains)
) {
	assert(next_chains_);
	auto &next_chains = *next_chains_;
	next_chains.clear();

	//For now, just test smoothed geodesic code to see if we can make some nicely spaced possible next lines.

	//based on example code in libgeodesic README:
	
	static hmContext &context = []() -> hmContext & {
		static hmContext context;
		hmContextInitialize( &context );
		return context;
	}();
	(void)context;

	hmTriMesh surface;
	hmTriDistance distance;

	hmTriMeshInitialize( &surface );
	hmTriDistanceInitialize( &distance );

	//translate model into libgeodesic's expected format:

	std::vector< double > verts;
	verts.reserve(model.vertices.size() * 3);
	for (auto const &v : model.vertices) {
		verts.emplace_back(v.x);
		verts.emplace_back(v.y);
		verts.emplace_back(v.z);
	}

	std::vector< size_t > faces;
	faces.reserve(model.triangles.size() * 3);
	for (auto const &f : model.triangles) {
		faces.emplace_back(f.x);
		faces.emplace_back(f.y);
		faces.emplace_back(f.z);
	}

	hmTriMeshReferenceData( &surface, verts.size() / 3, verts.data(), faces.size() / 3, faces.data() );
	distance.surface = &surface;

	hmTriDistanceEstimateTime( &distance );
	//TODO: potentially smooth distance somewhat by increasing distance->time
	std::cout << "time was set at " << distance.time;
	distance.time *= 1.0;
	std::cout << " heuristically, we've modified it to " << distance.time;

	hmTriDistanceSetBoundaryConditions( &distance, 0.5 );

	hmTriDistanceBuild( &distance );

	hmClearArrayDouble( distance.isSource.values, distance.surface->nVertices, 0.0 );

	//awkward -- should probably embed chains before doing this. Will work for testing for now, though:
	for (auto const &chain : active_chains) {
		for (auto const &ev : chain) {
			for (uint32_t i = 0; i < 3; ++i) {
				if (ev.simplex[i] != -1U) {
					distance.isSource.values[ev.simplex[i]] = 1.0;
				}
			}
		}
	}

	hmTriDistanceUpdate( &distance );

	//DEBUG: just show several distances:
	for (int step = 1; step < 10; ++step) {

		float level = step * (2.0f * parameters.stitch_height_mm / parameters.model_units_mm);

		std::vector< float > values;
		values.assign(distance.distance.values, distance.distance.values + distance.distance.nRows);
		assert(values.size() == model.vertices.size());

		std::vector< std::vector< EmbeddedVertex > > temp_chains;

		ak::extract_level_chains(model, values, level, &temp_chains);
		
		next_chains.insert(next_chains.end(), temp_chains.begin(), temp_chains.end());
	}

	hmTriDistanceDestroy( &distance );
	hmTriMeshDestroy( &surface );

}

