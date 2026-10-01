# Development image: the conda environment from environment.yml, with ioNERDSS
# installed in editable mode from /app, serving Jupyter on port 8888. README.md
# ("Docker development environment") has the build and run commands.

# Miniforge's conda is configured with conda-forge alone, the only channel
# environment.yml asks for, so the solve never reaches Anaconda's `defaults`
# channels. continuumio/miniconda3, used before, is deprecated.
FROM condaforge/miniforge3:latest

WORKDIR /app

# The source has to be in place before the environment is created:
# environment.yml pip-installs `-e .[test,jupyter]`, and conda runs pip from the
# directory holding the environment file, so `.` is /app and must already
# contain pyproject.toml. .dockerignore keeps .git and build output out.
COPY . .

# The name is set here, next to the PATH that depends on it, rather than taken
# from environment.yml.
RUN conda env create --name ionerdss-dev --file environment.yml && \
    conda clean --all --yes && \
    rm -rf /root/.cache/pip

# Commands run by the CMD, `docker run` and `docker exec` use the environment.
ENV PATH=/opt/conda/envs/ionerdss-dev/bin:$PATH

# Interactive shells (`docker exec -it ... bash`, Jupyter terminals) source
# ~/.bashrc, where the base image activates conda's base environment, which
# would put base's python ahead of this environment's.
RUN echo "conda activate ionerdss-dev" >> ~/.bashrc

# Fail the build, not the first `docker run`, if python or jupyter do not come
# from the environment. -I keeps /app off sys.path, so ioNERDSS has to be found
# through its editable install; the two attributes import its dependencies.
RUN python -I -c "import ionerdss; ionerdss.build_system_from_pdb; ionerdss.Analyzer" && \
    python -I -c "import importlib.metadata as m; print('ioNERDSS', m.version('ioNERDSS'))" && \
    jupyter notebook --version

EXPOSE 8888

CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
