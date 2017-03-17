package Webserver.Servlet;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;

/**
 * Created by antoine on 17/03/17.
 */
public class DownloadInputFileServlet extends HttpServlet {

    private String fileName;
    private String basePath;

    public DownloadInputFileServlet (String fileName){
        this.fileName = fileName;
        basePath = "/home/antoine/Documents/OSUR/Docker/hs1d/test_case/matlab/custom_input/";
    }

    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        response.setContentType("text/html");
        response.setStatus(HttpServletResponse.SC_OK);
        // TODO
    }

    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        response.setContentType("text/html"); // not sure
        response.setStatus(HttpServletResponse.SC_OK);
        // TODO
    }

}
