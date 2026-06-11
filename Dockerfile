# Utilisation d'une base Ubuntu récente
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    python3-pip \
    cargo \
    gcc \
    g++ \
    libomp-dev \
    bcalm \
    python3-numpy \
    python3-pysam \
    fastp \
    seqtk \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# 2. On copie TOUT le dossier de votre projet (script + sources C++/Rust)
COPY . /app/

# 3. Compilation du code Rust et C++ à la création de l'image Docker
# (Je remplace ${BIN_DIR} par le chemin relatif ./bin)
RUN cargo build --release --manifest-path ./bin/gene_finder_de_novo/Cargo.toml && \
    cp ./bin/gene_finder_de_novo/target/release/gene_finder_de_novo ./bin/gene_finder_de_novo.exe && \
    \
    cargo build --release --manifest-path ./bin/filtering_low_ab_percent/Cargo.toml && \
    cp ./bin/filtering_low_ab_percent/target/release/filtering_low_ab_percent ./bin/filtering_low_ab_percent.exe && \
    \
    cargo build --release --manifest-path ./bin/at_compressor/Cargo.toml && \
    cp ./bin/at_compressor/target/release/at_compressor ./bin/homopolymorphic_compression.exe && \
    \
    g++ -O3 -g -fopenmp ./bin/graph.cpp ./bin/ponderation.cpp -o ./bin/graph.exe && \
    g++ -O3 -g ./bin/graph.cpp ./bin/agglo.cpp -o ./bin/agglo.exe


# Rendre le script exécutable
RUN chmod +x /app/ET-core.sh

# Définition du script comme point d'entrée
ENTRYPOINT ["/app/ET-core.sh"]
