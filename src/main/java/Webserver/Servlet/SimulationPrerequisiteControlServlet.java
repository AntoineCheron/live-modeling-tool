package Webserver.Servlet;

import com.google.gson.Gson;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * Servlet that provides the list of element needed in order to run a simulation
 * Created by antoine on 20/03/17.
 */
public class SimulationPrerequisiteControlServlet extends HttpServlet {

    // Attributes
    private ArrayList<String> inputFileList;
    private ArrayList<String> missingFiles;

    // Constructor
    public SimulationPrerequisiteControlServlet() {
        inputFileList = new ArrayList<String>();
        inputFileList.add("geologic.input");
        inputFileList.add("morphologic.input");
        inputFileList.add("hydrologic.input");
        inputFileList.add("key_spatialized_param.param");
    }

    // HttpServlet methods

    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        // First, generate the list of missing files
        generateMissingInputList();

        // Configure response header
        response.setContentType("application/json");
        response.setStatus(HttpServletResponse.SC_OK);

        // Put the information into the response object
        PrintWriter out = response.getWriter();
        // Add the list into the out object after converting it into a json
        out.println(new Gson().toJson(missingFiles));
    }

    // Custom methods

    // Not working : TODO resolve it
    private void generateMissingInputList () throws NullPointerException {
        missingFiles = (ArrayList<String>) inputFileList.clone();

        // Path containing the files that will be used by the simulation
        String filePath = Utils.Path.getSimulationInputFilePath();

        // Retrieve the list of files contained into the folder
        File dir = new File(filePath);
        File[] dir_contents = dir.listFiles();

        /* Check whether each file is in the destination folder or not and if it finds a needed file, remove it
        from the list of missing files */
        if(dir_contents != null){
            for (File dir_content : dir_contents) {
                if (missingFiles.contains(dir_content.getName())) {
                    missingFiles.remove(dir_content.getName());
                }
            }
        }

    }
}
