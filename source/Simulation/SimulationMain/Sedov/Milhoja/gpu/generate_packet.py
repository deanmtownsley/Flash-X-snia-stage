import milhoja
from milhoja import DataPacketGenerator

if __name__ == "__main__":

    sizes = {
        "real": 8,
        "int": 4,
        "unsigned int": 4,
        "std::size_t": 8,
        "IntVect": 8,
        "RealVect": 16,
        "bool": 1
    }

    destination = "./"
    overwrite = True
    tf_spec = milhoja.TaskFunction.from_milhoja_json("gpu_tf_hydro_3D.json")
    logger = milhoja.BasicLogger(milhoja.LOG_LEVEL_BASIC)    

    generator = DataPacketGenerator(tf_spec, 4, logger, sizes)
    generator.generate_templates(destination, overwrite)
    generator.generate_source_code(destination, overwrite)
    generator.generate_header_code(destination, overwrite) 
