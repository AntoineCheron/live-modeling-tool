package Webserver;

import Webserver.Servlet.*;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.ServerConnector;
import org.eclipse.jetty.server.handler.ContextHandler;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.handler.ResourceHandler;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;

import javax.servlet.MultipartConfigElement;
import javax.servlet.http.HttpServlet;

/**
 *  Webserver set and configure the server
 *
 * Created by antoine on 17/03/17.
 */
public class Webserver {

    private static String inputFileBasePath = Utils.Path.getSimulationInputFilePath();
    public static void main (String[] args) throws Exception {
        // Define the server
        Server server = new Server();

        // Add http connector to server
        ServerConnector http = new ServerConnector(server);
        http.setHost("localhost");
        http.setPort(8080);
        http.setIdleTimeout(30000);
        server.addConnector(http);

        // Create a new RessourceHandler that serves the static content (html, css and js files) and configures it
        ResourceHandler resourceHandler = new ResourceHandler();

        ContextHandler contextHandler1 = new ContextHandler();
        contextHandler1.setContextPath("/");
        contextHandler1.setResourceBase("src_web/");
        contextHandler1.setWelcomeFiles(new String[]{"index.html"});
        contextHandler1.setHandler(resourceHandler);

        // Create a servlet context handler
        ServletContextHandler servlets = new ServletContextHandler(ServletContextHandler.SESSIONS);
        servlets.setContextPath("/API");

        // Configure all the API services based on servlets
        servlets.addServlet(getServletHolderForUpload(new UploadInputFileServlet("geologic.input"), inputFileBasePath), "/uploadGeolInput/*");
        servlets.addServlet(getServletHolderForUpload(new UploadInputFileServlet("hydrologic.input"), inputFileBasePath), "/uploadHydroInput/*");
        servlets.addServlet(getServletHolderForUpload(new UploadInputFileServlet("morphologic.input"), inputFileBasePath), "/uploadMorphoInput/*");
        servlets.addServlet(getServletHolderForUpload(new UploadInputFileServlet("key_spatialized_param.param"), inputFileBasePath), "/uploadParamInput/*");
        servlets.addServlet(new ServletHolder(new SimulationServlet()), "/simulate/*");
        servlets.addServlet(new ServletHolder(new SimulationPrerequisiteControlServlet()), "/missingElementToSimulate/*");
        servlets.addServlet(new ServletHolder(new SimulationResultsFormatServlet()), "/resultsFormat/*");
        servlets.addServlet(new ServletHolder(new SimulationResultsServlet()), "/result/*");


        // Add the contextHandler and the servlet context handler to the server
        HandlerList handlers = new HandlerList();
        handlers.setHandlers(new Handler[]{ contextHandler1, servlets});
        server.setHandler(handlers);

        // Start the server (because it is configured :) )
        // Server.join is used to make the server join the current thread
        server.start();
        server.join();
    }

    /* ================================================================================================================
     ==================================================================================================================
                                                    Custom methods
    ==================================================================================================================
    =================================================================================================================*/

    // Return a ServletHolder that accept file uploading for its servlet
    private static ServletHolder getServletHolderForUpload (HttpServlet servlet, String basePath) {
        // In order to create a servletHolder that accepts file upload, we need to :

        // Create a new servletHolder based on the servlet we want to use
        ServletHolder sHolder = new ServletHolder(servlet);

        // Configure the servlet holder to accept file uploading (using MultipartConfig)
        sHolder.getRegistration().setMultipartConfig(new MultipartConfigElement(inputFileBasePath));

        // Return it
        return sHolder;
    }
}
