package fr.irisa.diverse.Webserver.Servlet;

import com.google.gson.Gson;
import fr.irisa.diverse.Utils.Path;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Servlet that serves the list of files available from the simulation_results folder
 *
 * Created by antoine on 21/03/17.
 */
public class SimulationResultsFormatServlet extends HttpServlet {

    // Attributes
    private String pathToResults = Path.getSimulationResultsPath();

    // Constructor
    public SimulationResultsFormatServlet () {
        // Don't do anything for now
    }

    // HttpServlet method
    protected void doGet (HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        // Retrieve the list of results
        String results = getListOfPossibleResults();

        // Set the response header
        response.setStatus(HttpServletResponse.SC_OK);

        // Add the result into the response
        response.getWriter().println(results);

        // FINISHED :D
    }

    // Custom methods
    private String getListOfPossibleResults () throws NullPointerException{
        // Retrieve the list of files contained into the simulation_results folder
        File dir = new File(pathToResults);
        File[] dir_contents = dir.listFiles();

        // Put the name of each file into an ArrayList
        ArrayList<String> list = new ArrayList<String>();

        for (File dir_content : dir_contents != null ? dir_contents : new File[0]) {
            list.add(dir_content.getName());
        }

        // Convert the ArrayList into a JSON and Return this list
        return new Gson().toJson(list);
    }
}
