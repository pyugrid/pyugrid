from pyugrid import UGrid
import test_util

# list of dictionaries, each dictionary defines a grid. Use empty lists
# if not specifying a variable. i.e. 'edges': []
test_datasets = [
	{
		'nodes': [
			[0.,0.],
        	[2.,0.],
        	[1.,2.],
        	[3.,2.],
        	],
    	'faces': [
    		[0, 1, 2],
    		[1, 3, 2],
    		],
    	'edges': [
    		[0,1],
    		[1,3],
    		[3,2],
    		[2,0],
    		],
    },
]


def test_create_ugrid():
	for dataset in test_datasets:
		grid = UGrid(nodes=dataset['nodes'], faces=dataset['faces'], edges=dataset['edges'])
		assert grid.nodes.tolist() == dataset['nodes']
		assert grid.faces.tolist() == dataset['faces']
		assert grid.edges.tolist() == dataset['edges']