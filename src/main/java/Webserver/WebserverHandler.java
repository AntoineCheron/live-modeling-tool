package Webserver;

import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.handler.AbstractHandler;

import javax.servlet.ServletException;
import javax.servlet.ServletOutputStream;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * The handler will handle each request and return the right answer
 *
 * Created by antoine on 17/03/17.
 */
public class WebserverHandler extends AbstractHandler {

    // Attributes
    private String greeting;
    private String body;

    // Constructor
    WebserverHandler() {
        greeting = "Hello World from WebserverHandler";
        body = "Body";
    }

    // Custom method
    public void handle (String target, Request baseRequest, HttpServletRequest request,
                        HttpServletResponse response) throws IOException, ServletException {
        //Log some things to test
        System.out.println(baseRequest.getPathInfo());

        response.setContentType("text/html; charset=utf-8");
        response.setStatus(HttpServletResponse.SC_OK);

        PrintWriter out = response.getWriter();

        out.println("<h1>" + greeting + "</h1>");

        baseRequest.setHandled(true);
    }
}
