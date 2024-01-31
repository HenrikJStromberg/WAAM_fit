def test_total():
    import WAAMEvaluator
    input_path = "examples/Testobjekt_aehnlich_paper.stp"
    output_path = "output/test_WAAMEvaluator"
    step_size = 0.5

    WAAMEvaluator.evaluateSpheres(input_path, output_path, step_size)

def test_meshing():
    import WAAMEvaluator
    input_path = "examples/Testobjekt_aehnlich_paper.stp"
    step_size = 0.5

    return WAAMEvaluator.getTriangulation(input_path, step_size)

def test_heuver_spheres(mesh):
    import rust_methods
    import numpy as np
    nc, inz, centers, normals, elementTags = mesh
    r_inner = rust_methods.get_sphere_radii(centers, -normals, elementTags.tolist()) # type: ignore
    r_inner = np.array(r_inner)
    r_outer = rust_methods.get_sphere_radii(centers, normals, elementTags.tolist()) # type: ignore
    r_outer = np.array(r_outer)

if __name__ == '__main__':
    import timeit
    total_evaluation_time = timeit.timeit("test_total()", setup="from __main__ import test_total", number=1)
    meshing_time = timeit.timeit("test_meshing()", setup="from __main__ import test_meshing", number=1)
    
    mesh = test_meshing()
    heuver_time = timeit.timeit("test_heuver_spheres(mesh)", globals=globals(), number=1)

    print("\nTotal Evaluation Time: {:.4}s".format(total_evaluation_time))
    print("Meshing: {:.4}s ~ {:.2%}".format(meshing_time, meshing_time/total_evaluation_time))
    print("Radii evaluation: {:.4}s ~ {:.2%}".format(heuver_time, heuver_time/total_evaluation_time))