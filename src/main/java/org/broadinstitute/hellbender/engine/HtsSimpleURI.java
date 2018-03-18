package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.*;
import java.util.HashMap;

/**
 * Default implementation for HtsURI.
 */
public class HtsSimpleURI implements HtsURI {

    private final String uriString;
    private URI uri;
    private Path cachedPath;

    // TODO: should we use the normalized uri string for hash/compare purposes, or the user-supplied string ?
    // TODO: special handling for authority/credentials, query params, fragment?

    public HtsSimpleURI(final String uriString) {
        Utils.nonNull(uriString);
        this.uriString = uriString;

        try {
            uri = new URI(uriString);
        } catch (URISyntaxException e) {
            final String errorMessage = String.format("%s must be a valid URI. %s", uriString, e.getMessage());
            throw new IllegalArgumentException(errorMessage);
        }
    }

    /**
     * Converts the URI to a {@link Path} object. If the filesystem cannot be found in the usual way, then attempt
     * to load the filesystem provider using the thread context classloader. This is needed when the filesystem
     * provider is loaded using a URL classloader (e.g. in spark-submit).
     *
     * Also makes an attempt to interpret the argument as a file name if it's not a URI.
     *
     * @return the resulting {@code Path}
     * @throws UserException if an I/O error occurs when creating the file system
     */
    @Override
    public Path toPath() {
        if (cachedPath != null) {
            return cachedPath;
        }

        try {
            Path tmpPath;
            if (getURI().getScheme() == null) {
                tmpPath = Paths.get(getURIString());
            } else {
                tmpPath = Paths.get(getURI());
            }
            setCachedPath(tmpPath);
        } catch (FileSystemNotFoundException e) {
            // TODO: IOUtils creates a new FileSystem on FileSystemNotFoundException, but it results in
            // TODO: unknown providers being accepted
            try {
                ClassLoader cl = Thread.currentThread().getContextClassLoader();
                if (cl == null) {
                    throw e;
                }
                setCachedPath(FileSystems.newFileSystem(getURI(), new HashMap<>(), cl).provider().getPath(getURI()));
            } catch (ProviderNotFoundException pe) {
                // not a valid URI. Caller probably just gave us a file name.
                setCachedPath(Paths.get(getURIString()));
            } catch (IOException io) {
                setCachedPath(null);
                throw new UserException(getURIString() + " is not a supported path", io);
            }
        }
        return cachedPath;
    }

    @Override
    public URI getURI() {
        return uri;
    }

    @Override
    public String getURIString() {
        return uriString;
    }

    @Override
    public boolean isPath() {
        try {
            return getPath() != null || toPath() != null;
        } catch (ProviderNotFoundException| FileSystemNotFoundException | IllegalArgumentException | UserException e) {
            return false;
        }
    }

    @Override
    public String getToPathFailureReason() {
        try {
            toPath();
            return "Success";
        } catch (ProviderNotFoundException e) {
            return String.format("ProviderNotFoundException: %s", e.getMessage());
        } catch (FileSystemNotFoundException e) {
            return String.format("FileSystemNotFoundException: %s", e.getMessage());
        } catch (IllegalArgumentException e) {
            return String.format("IllegalArgumentException: %s", e.getMessage());
        } catch (UserException e) {
            return String.format("UserException: %s", e.getMessage());
        }
    }

    // get the path associated with this URI
    protected Path getPath() { return cachedPath; }

    protected void setCachedPath(Path path) {
        this.cachedPath = path;
    }

    @Override
    public String toString() {
        return uriString;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof HtsSimpleURI)) return false;

        HtsSimpleURI that = (HtsSimpleURI) o;

        if (!getURIString().equals(that.getURIString())) return false;
        if (!getURI().equals(that.getURI())) return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = getURIString().hashCode();
        result = 31 * result + getURI().hashCode();
        return result;
    }
}
