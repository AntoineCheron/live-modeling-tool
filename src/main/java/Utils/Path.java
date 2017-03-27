package Utils;

/**
 * This file is a util, used to avoid hard coding path to simulation files
 * Created by antoine on 21/03/17.
 */
public class Path {

    public static String getHs1dPath () {
        return System.getProperty("user.dir") + "/SimulationUtils/";
    }

    public static String getSimulationInputFilePath () {
        return System.getProperty("user.dir") + "/SimulationUtils/hs1d/test_case/matlab/custom_input/";
    }

    public static String getSimulationResultsPath () {
        return System.getProperty("user.dir") + "/SimulationUtils/hs1d/src_python/simulation_results/";
    }
}
