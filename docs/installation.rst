Installation
============

``snps`` is `available <https://pypi.org/project/snps/>`_ on the
`Python Package Index <https://pypi.org>`_. Install ``snps`` (and its required
Python dependencies) via ``pip``::

    $ pip install snps

Installation and Usage on a Raspberry Pi
----------------------------------------
The instructions below provide the steps to install ``snps`` on a
`Raspberry Pi <https://www.raspberrypi.org>`_ (tested with
"`Raspberry Pi OS <https://www.raspberrypi.org/downloads/raspberry-pi-os/>`_ (32-bit) Lite",
release date 2020-08-20). For more details about Python on the Raspberry Pi, see
`here <https://www.raspberrypi.org/documentation/linux/software/python.md>`_.

.. note:: Text after a prompt (e.g., ``$``) is the command to type at the command line. The
          instructions assume a fresh install of Raspberry Pi OS and that after logging in as
          the ``pi`` user, the current working directory is ``/home/pi``.

1. Install ``pip`` for Python 3::

    pi@raspberrypi:~ $ sudo apt install python3-pip

   Press "y" followed by "enter" to continue. This enables us to install packages from the
   Python Package Index.

2. Install the ``venv`` module::

    pi@raspberrypi:~ $ sudo apt install python3-venv

   Press "y" followed by "enter" to continue. This enables us to create a
   `virtual environment <https://docs.python.org/3/library/venv.html>`_ to isolate the ``snps``
   installation from other system Python packages.

3. `Install ATLAS <https://github.com/Kitt-AI/snowboy/issues/262#issuecomment-324997127>`_::

    pi@raspberrypi:~ $ sudo apt install libatlas-base-dev

   Press "y" followed by "enter" to continue. This is required for `NumPy <https://numpy.org>`_, a
   dependency of ``snps``.

4. Create a directory for ``snps`` and change working directory::

    pi@raspberrypi:~ $ mkdir snps
    pi@raspberrypi:~ $ cd snps

5. Create a virtual environment for ``snps``::

    pi@raspberrypi:~/snps $ python3 -m venv .venv

   The virtual environment is located at ``/home/pi/snps/.venv``.

6. Activate the virtual environment::

    pi@raspberrypi:~/snps $ source .venv/bin/activate

   Now when you invoke Python or ``pip``, the virtual environment's version will be used (as
   indicated by the ``(.venv)`` before the prompt). This can be verified as follows::

    (.venv) pi@raspberrypi:~/snps $ which python
    /home/pi/snps/.venv/bin/python

7. Install ``snps``::

    (.venv) pi@raspberrypi:~/snps $ pip install snps

8. Start Python::

    (.venv) pi@raspberrypi:~/snps $ python
    Python 3.7.3 (default, Jul 25 2020, 13:03:44)
    [GCC 8.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

9. Use ``snps``; examples shown in the README should now work.

10. At completion of usage, the virtual environment can be deactivated::

     (.venv) pi@raspberrypi:~/snps $ deactivate
     pi@raspberrypi:~/snps $

