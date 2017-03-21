package Webserver.Servlet;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.Part;
import java.io.IOException;

/**
 * Servlet used to let the user upload the files that are needed to run a simulation
 * Created by antoine on 17/03/17.
 */
public class UploadInputFileServlet extends HttpServlet {

    private String fileName;

    public UploadInputFileServlet(String fileName){
        this.fileName = fileName;
    }

    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

        // Retrieve the file
        Part file = request.getPart("file");

        // Save the file onto the hard drive with the right name
        file.write(fileName);

        // Set the status to tell the client everything went well
        response.setStatus(HttpServletResponse.SC_OK);
    }

}
