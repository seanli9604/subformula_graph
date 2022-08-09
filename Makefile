build:
	DOCKER_BUILDKIT=1 docker build -t subformula_graph_app . 

run:
	docker run subformula_graph_app