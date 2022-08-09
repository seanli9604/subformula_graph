build:
	DOCKER_BUILDKIT=1 docker build -t subformula_graph_app . 

run:
	docker run -p 80:5000 subformula_graph_app 