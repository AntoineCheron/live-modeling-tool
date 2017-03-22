package Webserver.Servlet;

import com.google.gson.Gson;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.*;
import java.nio.file.Files;
import java.util.List;

/**
 * This servlet serves the content of the simulation's result files as Json object.
 *
 * Created by antoine on 21/03/17.
 */
public class SimulationResultsServlet extends HttpServlet {

    // Attributes
    private String pathToResults = Utils.Path.getSimulationResultsPath();

    // Constructor
    public SimulationResultsServlet () {
        // Don't do anything for now
    }

    // HttpServlet method
    protected void doGet (HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException, NullPointerException {
        // Retrieve the file name (= resultName) from the request
        String name = request.getParameter("name");

        // Verify that a name has been sent
        if(name == null) {
            // Tell the client that the request hasn't been correctly formatted
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        }
        else {
            // Retrieve the file containing the requested result
            File results_dir = new File(pathToResults);
            File[] dir_content = results_dir.listFiles();

            // Verify that this file exist. If not : return a BAD_REQUEST response
            boolean fileInFolder = false;
            for (File aDir_content : dir_content) {
                if (aDir_content.toString().equals(pathToResults + name)) {
                    fileInFolder = true;
                    break;
                }
            }
            // Check the boolean
            if(!fileInFolder) {
                response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
                return;
            }

            // If the file has been found into the folder, we can continue
            File resultFile = new File(pathToResults + name);

            // Parse the file to get the correctly formatted json
            String result = convertToJson(resultFile);

            // File successfully read and formatted, time to build the response
            response.setStatus(HttpServletResponse.SC_OK);
            response.setContentType("application/json");
            response.getWriter().println(result);
        }
    }

    // Custom methods
    private String[][] convertFileIntoStringArray (File file) throws IOException {
        List<String> lines = Files.readAllLines(file.toPath());

        // The result will be an array of string with the same amount of line as the original file
        String[] tempResult = new String[lines.size()];

        // Transfer content of lines into result
        for(int i=0; i<lines.size(); i++){
            tempResult[i] = lines.get(i);
        }

        // Now, parse each line to separate the columns in a result[line][column]

        // First : need to know how much columns there are in the file. The information is in firstLineWithColumns.length
        String[] firstLineWithColumns = tempResult[0].split("\t");

        String[][] result = new String[tempResult.length][firstLineWithColumns.length];

        // Split each line into the proper number of columns, splitting the string for each space
        for (int line=0; line<tempResult.length; line++) {
            result[line] = tempResult[line].split("\t");
        }

        // Return the result, FINISH :D
        return result;
    }

    private String convertToJson (File file) throws IOException{

        String[][] terribleString = convertFileIntoStringArray(file);

        if(terribleString[0].length != 1) return new Gson().toJson(terribleString);
        // In case the terribleString array contains only one column, put it into a String[] instead of a String[][]
        else {
            String[] result = new String[terribleString.length];
            for(int i=0; i<terribleString.length; i++){
                result[i] = terribleString[i][0];
            }
            return new Gson().toJson(result);
        }


    }

}
