import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate, Link } from 'react-router-dom';

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from 'uuid'

import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-enterprise';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';
import './ag-theme-acmecorp.css';

import Select from 'react-select';

import { variables, AG_GRID_LOCALE_RU } from './Variables.js';


var ErrorMessage = 0
const per_topics = ['Поиск в pubmed', 'Тематический анализ', 'Поиск в векторном представлении']


export class AdminPanel extends Component {

  constructor(props) {
    super(props);

    this.gridRef = createRef();
    this.searchgridRef = createRef();
    this.analisegridRef = createRef();
    this.state = {
      token: variables.token,
      loading: false,

      useUsers: true,
      useSearch: false,
      useAnalise: false,

      permissions: [],
      users: null,
      DetailUser: null,
      usersInfo: [
        { field: 'email', filter: 'agTextColumnFilter', enableRowGroup: true, minWidth: 100, width: 300, resizable: true},
        { field: 'allow_status', filter: 'agNumberColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
        { field: 'is_admin', sortable: true, enableRowGroup: true, resizable: true},
        { field: 'last_login', sortable: true, enableRowGroup: true, resizable: true}
      ],
      search_queries: null,
      DetailSearch: null,
      searchInfo: [
        { field: 'query', filter: 'agTextColumnFilter', enableRowGroup: true, minWidth: 100, width: 300, resizable: true},
        { field: 'full_query', filter: 'agTextColumnFilter', enableRowGroup: true, minWidth: 100, width: 300, resizable: true},
        { field: 'translation_stack', filter: 'agTextColumnFilter', enableRowGroup: true, minWidth: 100, width: 300, resizable: true},
        { field: 'status', filter: 'agNumberColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
        { field: 'user', filter: 'agTextColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
        { field: 'start_date', filter: 'agTextColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
        { field: 'work_time', filter: 'agNumberColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
      ],
      analise_queries: null,
      DetailAnalise: null,
      analiseInfo: [
        { field: 'type_analise', filter: 'agTextColumnFilter', enableRowGroup: true, minWidth: 100, width: 300, resizable: true},
        { field: 'status', filter: 'agTextColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
        { field: 'user', filter: 'agTextColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
        { field: 'start_date', filter: 'agTextColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
        { field: 'work_time', filter: 'agNumberColumnFilter', sortable: true, enableRowGroup: true, resizable: true},
      ],

      Email: '',
      Password: '',
      DualPassword: '',

      max_search: 9999,
      max_analise: 9999,
      max_ddi: 9999,
      changedUser: null,
      allow_type: {label: 'Все', value: 3},
      allow_types: [{label: 'Тематический анализ', value: 0}, {label: 'Факты EBM', value: 1}, {label: 'Поговорим', value: 2}, {label: 'Все', value: 3}]
    }
  }

   // Permissions
  getPermissions() {
    fetch(variables.API_URL + '/api/permissions',
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then(response => {
        console.log(response.status);
        ErrorMessage = response.status
        if (response.ok) {
          return response.json()
        } else {
          throw Error(response.status)
        }
      })
      .then(data => {
        console.log(data)
        this.setState({permissions: data.permissions})
      })
      .catch(error => {
        console.log(error)
      })
  }

  getInfo() {
    fetch(variables.API_URL + '/api/admin', {
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      }
    })
    .then((res) => {
        if (res.ok) { return res.json() }
        else { throw Error(res.statusText) }
      })
      .then((result) => {
        console.log(result.data)
        this.setState({
            users: result.users, DetailUser: result.users[0],
            search_queries: result.search, DetailSearch: result.search[0],
            analise_queries: result.analise, DetailAnalise: result.analise[0]
        })
      })
      .catch((error) => {
        console.log(error)
      })
  }

  componentDidMount() {
    this.getInfo();
    this.getPermissions();
    console.log('start');
  }

  onSelectionChanged = () => {
    const selectedRows = this.gridRef.current.api.getSelectedRows();
    this.setState({ DetailUser: (selectedRows.length === 1 ? selectedRows[0] : null) })
  }

  onSelectionSearchChanged = () => {
    const selectedRows = this.searchgridRef.current.api.getSelectedRows();
    this.setState({ DetailSearch: (selectedRows.length === 1 ? selectedRows[0] : null) })
  }

  onSelectionAnaliseChanged = () => {
    const selectedRows = this.analisegridRef.current.api.getSelectedRows();
    this.setState({ DetailAnalise: (selectedRows.length === 1 ? selectedRows[0] : null) })
  }

  setUserModal() {
    this.setState({ Email: '', Password: '', DualPassword: '', message: '', messageStatus: 0 })
  }

  setUpdateModal(user) {
    this.setState({
        max_search: 9999,
        max_analise: 9999,
        max_ddi: 9999,
        allow_type: {label: 'Все', value: 2},
        changedUser: user,
        message: '',
        messageStatus: 0
    })
  }

  createUserClick() {
    if (this.state.Email === '') {
        alert('Введите email!')
        return
    }
    if (this.state.Password === '' || this.state.DualPassword === '') {
        alert('Введите пароль!')
        return
    }
    if (this.state.Password !== this.state.DualPassword) {
        alert('Пароли не совпадают!')
        return
    }
    fetch(variables.API_URL + '/accounts/create', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        email: this.state.Email,
        password: this.state.Password,
      })
    })
      .then((res) => {
        if (res.ok) { return res.json() }
        else { throw Error(res.statusText) }
      })
      .then((result) => {
        this.setState({ message: 'Пользователь создан', messageStatus: 201, users: result.users, DetailUser: result.users[0]})
      })
      .catch((error) => {
        this.setState({ message: 'Ошибка при создании', messageStatus: 500})
      })
  }

  updateUserClick() {
    fetch(variables.API_URL + '/accounts/create', {
      method: 'PUT',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        user: this.state.changedUser.email,
        allow_type: this.state.allow_type.value,
        0: this.state.max_search,
        1: this.state.max_analise,
        2: this.state.max_ddi
      })
    })
      .then((res) => {
        if (res.ok) { return res.json() }
        else { throw Error(res.statusText) }
      })
      .then((result) => {
        this.setState({ message: 'Пользователь обновлен', messageStatus: 201, DetailUser: result.user })
        this.gridRef.current.api.forEachNode((rowNode) => {
            if (rowNode.data.email !== result.user.email) {
              return;
            }
            console.log('some')
            // arbitrarily update some data
            const updated = result.user;
            // directly update data in rowNode
            rowNode.updateData(updated);
          })
      })
      .catch((error) => {
        this.setState({ message: 'Ошибка при обновлении', messageStatus: 500})
      })
  }

  render() {
    const {
      token,
      loading,
      permissions,

      users,
      usersInfo,
      DetailUser,

      search_queries,
      searchInfo,
      DetailSearch,

      analise_queries,
      analiseInfo,
      DetailAnalise,

      useUsers,
      useSearch,
      useAnalise,

      Email,
      Password,
      DualPassword,

      message,
      messageStatus,

      max_analise,
      max_search,
      max_ddi,
      allow_type,
      allow_types,

    } = this.state;

    if (!token) {
      return <Navigate push to="/login" />
    } else {
      return (
        <>
          <header>
            <nav class="bg-white border-gray-200 px-4 lg:px-6 py-2.5">
              <div class="flex flex-wrap justify-between items-center">
                <div class="flex justify-start items-center">
                  <a href="" class="flex mr-4">
                    <img src="https://flowbite.s3.amazonaws.com/logo.svg" class="mr-3 h-8" alt="FlowBite Logo" />
                    <span class="self-center text-2xl font-semibold whitespace-nowrap">EBM Sechenov DataMed.AI</span>
                  </a>
                    <ul class="flex font-medium flex-row space-x-8">
                      <Link to="/tematic_review">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Тематический анализ</a>
                        </li>
                      </Link>
                      <Link to="/ddi_review">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Факты для EBM</a>
                        </li>
                      </Link>
                      <Link to="/admin">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 bg-blue-700 rounded md:bg-transparent md:text-blue-700 md:p-0" aria-current="page">Админ панель</a>
                        </li>
                      </Link>
                    </ul>
                </div>
                <div class="flex items-center lg:order-3">
                  <div class="flex-shrink-0 dropdown">
                    <a href="#" class="d-block link-body-emphasis text-decoration-none dropdown-toggle" data-bs-toggle="dropdown" aria-expanded="false">
                      <img src="https://github.com/mdo.png" alt="mdo" width="32" height="32" class="rounded-circle" />
                    </a>
                    <ul class="dropdown-menu text-small shadow">
                      {permissions?.map(per =>
                        <li><a class="dropdown-item" href="#">{per.topic} {per.all_records? `${per.all_records}`: 'безлимитно'}</a></li>
                      )}
                    </ul>
                  </div>
                </div>
                <div class="flex items-center lg:order-2">
                  <button type="button" class="hidden sm:inline-flex items-center justify-center text-white bg-primary-700 hover:bg-primary-800 focus:ring-4 focus:ring-primary-300 font-medium rounded-lg text-xs px-3 py-1.5 mr-2 dark:bg-primary-600 dark:hover:bg-primary-700 focus:outline-none dark:focus:ring-primary-800"><svg aria-hidden="true" class="mr-1 -ml-1 w-5 h-5" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M10 5a1 1 0 011 1v3h3a1 1 0 110 2h-3v3a1 1 0 11-2 0v-3H6a1 1 0 110-2h3V6a1 1 0 011-1z" clip-rule="evenodd"></path></svg> Действие</button>
                </div>
              </div>
            </nav>
            <nav class="bg-white border-gray-200 px-6">
              <div class="w-full">
                <div class="flex justify-between items-center">
                  <button id="toggleSidebar" aria-expanded="true" aria-controls="sidebar" class="hidden p-2 mr-3 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700" data-bs-toggle="collapse" data-bs-target="#sidebar" aria-label="Toggle navigation">
                    <svg class="w-6 h-6" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clip-rule="evenodd"></path></svg>
                  </button>
                  <button id="toggleSidebar" aria-expanded="true" aria-controls="sidebar2" class="order-last hidden p-2 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700" data-bs-toggle="collapse" data-bs-target="#sidebar2" aria-label="Toggle navigation">
                    <svg class="w-6 h-6 rotate-180" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clip-rule="evenodd"></path></svg>
                  </button>
                </div>
              </div>
            </nav>
          </header >
          <main>
            <div>
              <div className="container-fluid">
                <div className="row align-items-stretch b-height">
                  <aside id="sidebar" className="h-screen col-md-2 my-3 bg-white collapse show width border rounded-3 g-0">
                    <div className="accordion accordion-flush" id="accordionFlushExample">
                        <div class="ml-5 grow items-center justify-between hidden w-full md:flex md:w-auto md:order-1" id="navbar-sticky">
                            <ul class="nav nav-pills" id="myTab" role="tablist">
                              <li class="nav-item mr-2" role="presentation">
                                <button class="nav-link inline-block px-4 py-2 rounded-lg hover:text-gray-900 hover:bg-gray-100 active" id="home-tab" data-bs-toggle="tab" data-bs-target="#home" type="button" role="tab" aria-controls="home" aria-selected={useUsers} onClick={() => this.setState({useUsers: true, useSearch: false, useAnalise: false})}>Пользователи</button>
                              </li>
                              <li class="nav-item mr-2" role="presentation">
                                <button class="nav-link inline-block px-4 py-2 rounded-lg hover:text-gray-900 hover:bg-gray-100" id="profile-tab" data-bs-toggle="tab" data-bs-target="#profile" type="button" role="tab" aria-controls="profile" aria-selected={useSearch} onClick={() => this.setState({useUsers: false, useSearch: true, useAnalise: false})}>Запросы на поиск</button>
                              </li>
                              <li class="nav-item mr-2" role="presentation">
                                <button class="nav-link inline-block px-4 py-2 rounded-lg hover:text-gray-900 hover:bg-gray-100" id="contact-tab" data-bs-toggle="tab" data-bs-target="#contact" type="button" role="tab" aria-controls="contact" aria-selected={useAnalise} onClick={() => this.setState({useUsers: false, useSearch: false, useAnalise: true})}>Запросы на анализ</button>
                              </li>
                            </ul>
                        </div>
                    </div>
                  </aside>
                  <section class="col p-3 m-3 border rounded-3 bg-white overflow-auto">
                    <div class="accordion accordion-flush" id="accordion">
                      <div class="accordion-item">
                        <div id="flush-collapseSeven" class="collapse multi-collapse" aria-labelledby="flush-headingSeven" data-bs-target="#accordionFlushExample">

                        </div>
                      </div>
                      <div>
                        </div>
                    </div>
                    <div>
                      <div class="bd-example">
                        <div class="tab-content" id="myTabContent">
                          <div class="tab-pane fade active show" id="home" role="tabpanel" aria-labelledby="home-tab">
                            <div class="container-fluid g-0">
                              <div className="ag-theme-alpine ag-theme-acmecorp" style={{ height: 700 }}>
                                <AgGridReact
                                  ref={this.gridRef}
                                  rowData={users}
                                  columnDefs={usersInfo}
                                  pagination={true}
                                  rowSelection={'single'}
                                  onSelectionChanged={this.onSelectionChanged}
                                  localeText={AG_GRID_LOCALE_RU}
                                  sideBar={{
                                    toolPanels: [
                                      {
                                        id: 'columns',
                                        labelDefault: 'Columns',
                                        labelKey: 'columns',
                                        iconKey: 'columns',
                                        toolPanel: 'agColumnsToolPanel',
                                        minWidth: 225,
                                        width: 225,
                                        maxWidth: 225,
                                      },
                                      {
                                        id: 'filters',
                                        labelDefault: 'Filters',
                                        labelKey: 'filters',
                                        iconKey: 'filter',
                                        toolPanel: 'agFiltersToolPanel',
                                        minWidth: 180,
                                        maxWidth: 400,
                                        width: 250,
                                      },
                                    ],
                                    position: 'left',

                                  }}
                                >
                                </AgGridReact>
                              </div>
                            </div>
                            <div class="flex justify-center my-6">
                                <button
                                    type="button"
                                    className="flex text-white bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:ring-blue-300 font-medium rounded-lg text-sm px-5 py-2.5 ml-2 dark:bg-blue-600 dark:hover:bg-blue-700 focus:outline-none dark:focus:ring-blue-800 float-end"
                                    data-bs-toggle="modal"
                                    data-bs-target="#userModal"
                                    onClick={() => this.setUserModal()}>
                                    Создать пользователя
                                </button>
                            </div>
                          </div>
                          <div class="tab-pane fade" id="profile" role="tabpanel" aria-labelledby="profile-tab">
                            <div class="container-fluid g-0">
                              <div className="ag-theme-alpine ag-theme-acmecorp" style={{ height: 700 }}>
                                <AgGridReact
                                  ref={this.searchgridRef}
                                  rowData={search_queries}
                                  columnDefs={searchInfo}
                                  pagination={true}
                                  rowSelection={'single'}
                                  onSelectionChanged={this.onSelectionSearchChanged}
                                  localeText={AG_GRID_LOCALE_RU}
                                  sideBar={{
                                    toolPanels: [
                                      {
                                        id: 'columns',
                                        labelDefault: 'Columns',
                                        labelKey: 'columns',
                                        iconKey: 'columns',
                                        toolPanel: 'agColumnsToolPanel',
                                        minWidth: 225,
                                        width: 225,
                                        maxWidth: 225,
                                      },
                                      {
                                        id: 'filters',
                                        labelDefault: 'Filters',
                                        labelKey: 'filters',
                                        iconKey: 'filter',
                                        toolPanel: 'agFiltersToolPanel',
                                        minWidth: 180,
                                        maxWidth: 400,
                                        width: 250,
                                      },
                                    ],
                                    position: 'left',

                                  }}
                                >
                                </AgGridReact>
                              </div>
                            </div>
                          </div>
                          <div class="tab-pane fade" id="contact" role="tabpanel" aria-labelledby="contact-tab">
                            <div class="container-fluid g-0">
                              <div className="ag-theme-alpine ag-theme-acmecorp" style={{ height: 700 }}>
                                <AgGridReact
                                  ref={this.analisegridRef}
                                  rowData={analise_queries}
                                  columnDefs={analiseInfo}
                                  pagination={true}
                                  rowSelection={'single'}
                                  onSelectionChanged={this.onSelectionAnaliseChanged}
                                  localeText={AG_GRID_LOCALE_RU}
                                  sideBar={{
                                    toolPanels: [
                                      {
                                        id: 'columns',
                                        labelDefault: 'Columns',
                                        labelKey: 'columns',
                                        iconKey: 'columns',
                                        toolPanel: 'agColumnsToolPanel',
                                        minWidth: 225,
                                        width: 225,
                                        maxWidth: 225,
                                      },
                                      {
                                        id: 'filters',
                                        labelDefault: 'Filters',
                                        labelKey: 'filters',
                                        iconKey: 'filter',
                                        toolPanel: 'agFiltersToolPanel',
                                        minWidth: 180,
                                        maxWidth: 400,
                                        width: 250,
                                      },
                                    ],
                                    position: 'left',
                                    filters: true
                                  }}
                                >
                                </AgGridReact>
                              </div>
                            </div>
                          </div>
                        </div>
                      </div>
                    </div>
                  </section>

                  <aside id="sidebar2" class="col-md-4 h-screen collapse show width col p-3 my-3 border rounded-3 bg-white">
                    {useUsers?
                    <>
                    <h3 class="pb-2 mb-3 border-bottom">Подробное описание Пользователя</h3>
                        <nav class="small" id="toc">
                            {DetailUser ?
                            <div class="card mb-3">
                              <div class="card-body">
                                <p class="card-text">---------------------------------- </p>
                                <p class="card-text">Email : {DetailUser.email}</p>
                                <p class="card-text">---------------------------------- </p>
                                <p class="card-text">Администратор : {DetailUser.is_admin? <input class="form-check-input" type="checkbox" value="" id="flexCheckDefault" checked /> :null}</p>
                                <p class="card-text">Доступ : {DetailUser.allow_status}</p>
                                <p class="card-text">---------------------------------- </p>
                                {DetailUser.permissions?.map(per =>
                                    <>
                                        <p class="card-text"><small class="text-success">{per.topic} : {per.all_records? `${per.all_records}`: 'безлимитно'} </small></p>
                                        <p class="card-text"><small class="text-success">Начало использования: {per.start_time}</small></p>
                                        <p class="card-text"><small class="text-success">---------------------------------- </small></p>
                                    </>
                                )}
                              </div>
                              <input
                                className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2"
                                value="Обновить"
                                type='submit'
                                data-bs-toggle="modal"
                                data-bs-target="#updateModal"
                                onClick={() => this.setUpdateModal(DetailUser)} />
                            </div>
                            : null}
                        </nav>
                    </>
                    : useSearch?
                    <>
                        <h3 class="pb-2 mb-3 border-bottom">Подробное описание Поискового запроса</h3>
                        <nav class="small" id="toc">
                            {DetailSearch ?
                            <div class="card mb-3">
                              <div class="card-body">
                                <p class="card-text">Query :  {DetailSearch.query} (Find {DetailSearch.count})</p>
                                <p class="card-text">---------------------------------- </p>
                                <p class="card-text">Full query :  </p>
                                <p class="card-text">{DetailSearch.full_query}</p>
                                <p class="card-text">---------------------------------- </p>
                                <p class="card-text">Translation stack :  {DetailSearch.translation_stack}</p>
                                <p class="card-text">---------------------------------- </p>
                                <p class="card-text"><small class="text-success">Дата начала запроса : {DetailSearch.start_date} </small></p>
                                <p class="card-text"><small class="text-success">Дата конца запроса : {DetailSearch.end_date}</small></p>
                                <p class="card-text"><small class="text-success">Время выполнения : {DetailSearch.work_time} сек.</small></p>
                                <p class="card-text"><small class="text-success">Пользователь : {DetailSearch.user}</small></p>
                                <p class="card-text"><small class="text-success">Статус запроса : {DetailSearch.status} </small></p>
                              </div>
                            </div>
                            : null}
                        </nav>
                    </>
                    : useAnalise?
                    <>
                        <h3 class="pb-2 mb-3 border-bottom">Подробное описание запроса на анализ</h3>
                        <nav class="small" id="toc">
                            {DetailAnalise ?
                            <div class="card mb-3">
                              <div class="card-body">
                                <p class="card-text">Type analise :  {DetailAnalise.type_analise}</p>
                                <p class="card-text">---------------------------------- </p>
                                <p class="card-text"><small class="text-success">Дата начала запроса : {DetailAnalise.start_date} </small></p>
                                <p class="card-text"><small class="text-success">Дата конца запроса : {DetailAnalise.end_date}</small></p>
                                <p class="card-text"><small class="text-success">Время выполнения : {DetailAnalise.work_time} сек.</small></p>
                                <p class="card-text"><small class="text-success">Пользователь : {DetailAnalise.user}</small></p>
                                <p class="card-text"><small class="text-success">Статус запроса : {DetailAnalise.status} </small></p>
                              </div>
                            </div>
                            : null}
                        </nav>
                    </>
                    :null}
                  </aside>

                </div>
              </div>
            </div>
          </main>
          <div className="modal fade" id="userModal" tabIndex="-1" aria-hidden="true">
            <div className="modal-dialog modal-lg modal-dialog-centered">
                <div className="modal-content bg-slate-50 rounded-lg drop-shadow-md dark:bg-gray-800">
                    <div className="modal-header">
                        <h2 className="modal-title font-semibold text-gray-900 dark:text-white">Создать пользователя</h2>
                        <button type="button" class="absolute top-3 right-2.5 text-gray-400 bg-transparent hover:bg-gray-200 hover:text-gray-900 rounded-lg text-sm p-1.5 ml-auto inline-flex items-center dark:hover:bg-gray-800 dark:hover:text-white" data-modal-hide="authentication-modal">
                            <svg aria-hidden="true" class="w-5 h-5" data-bs-dismiss="modal" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M4.293 4.293a1 1 0 011.414 0L10 8.586l4.293-4.293a1 1 0 111.414 1.414L11.414 10l4.293 4.293a1 1 0 01-1.414 1.414L10 11.414l-4.293 4.293a1 1 0 01-1.414-1.414L8.586 10 4.293 5.707a1 1 0 010-1.414z" clip-rule="evenodd"></path></svg>
                            <span class="sr-only">Закрыть диалог</span>
                        </button>
                    </div>

                    <div className="modal-body">
                        <div class="relative w-full mb-6 group">
                            <label for="text" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white">Email</label>
                            <input
                                type="email"
                                className="bg-gray-50 border border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5 dark:bg-gray-700 dark:border-gray-600 dark:placeholder-gray-400 dark:text-white dark:focus:ring-blue-500 dark:focus:border-blue-500"
                                value={Email}
                                onChange={(e) => this.setState({Email: e.target.value})}
                                placeholder="email"
                            />
                            <label for="text" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white">Password</label>
                            <input
                                type="password"
                                className="bg-gray-50 border border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5 dark:bg-gray-700 dark:border-gray-600 dark:placeholder-gray-400 dark:text-white dark:focus:ring-blue-500 dark:focus:border-blue-500"
                                value={Password}
                                onChange={(e) => this.setState({Password: e.target.value})}
                                placeholder="**********"
                            />
                            <label for="text" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white">Dublicate Password</label>
                            <input
                                type="password"
                                className="bg-gray-50 border border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5 dark:bg-gray-700 dark:border-gray-600 dark:placeholder-gray-400 dark:text-white dark:focus:ring-blue-500 dark:focus:border-blue-500"
                                value={DualPassword}
                                onChange={(e) => this.setState({DualPassword: e.target.value})}
                                placeholder="**********"
                            />
                            <br />
                            {message ?
                                messageStatus > 299 ?
                                  <p class="pb-2 mb-3 border-bottom" style={{ color: 'red' }}>{message}.</p>
                                  : messageStatus === 201 ?
                                    <p class="pb-2 mb-3 border-bottom" style={{ color: 'green' }}>{message}.</p>
                                  :null
                            : null}
                            <button type="button"
                                className="flex text-white bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:ring-blue-300 font-medium rounded-lg text-sm px-5 py-2.5 mr-2 dark:bg-blue-600 dark:hover:bg-blue-700 focus:outline-none dark:focus:ring-blue-800"
                                onClick={() => this.createUserClick()}
                            >Создать</button>
                        </div>
                    </div>

                </div>
            </div>
          </div>
          <div className="modal fade" id="updateModal" tabIndex="-1" aria-hidden="true">
            <div className="modal-dialog modal-lg modal-dialog-centered">
                <div className="modal-content bg-slate-50 rounded-lg drop-shadow-md dark:bg-gray-800">
                    <div className="modal-header">
                        <h2 className="modal-title font-semibold text-gray-900 dark:text-white">Обновить доступ пользователя</h2>
                        <button type="button" class="absolute top-3 right-2.5 text-gray-400 bg-transparent hover:bg-gray-200 hover:text-gray-900 rounded-lg text-sm p-1.5 ml-auto inline-flex items-center dark:hover:bg-gray-800 dark:hover:text-white" data-modal-hide="authentication-modal">
                            <svg aria-hidden="true" class="w-5 h-5" data-bs-dismiss="modal" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M4.293 4.293a1 1 0 011.414 0L10 8.586l4.293-4.293a1 1 0 111.414 1.414L11.414 10l4.293 4.293a1 1 0 01-1.414 1.414L10 11.414l-4.293 4.293a1 1 0 01-1.414-1.414L8.586 10 4.293 5.707a1 1 0 010-1.414z" clip-rule="evenodd"></path></svg>
                            <span class="sr-only">Закрыть диалог</span>
                        </button>
                    </div>

                    <div className="modal-body">
                        <div class="relative w-full mb-6 group">
                            <label for="text" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white">Доступно записей для поиска в PubMed</label>
                            <input
                                type="number"
                                id="replyNumber1"
                                min="0"
                                step="1"
                                value={max_search}
                                data-bind="value:replyNumber1"
                                onChange={(e) => this.setState({max_search: e.target.value})}
                            />
                            <label for="text" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white">Доступно записей для тематического анализа</label>
                            <input
                                type="number"
                                id="replyNumber2"
                                min="0"
                                step="1"
                                value={max_analise}
                                data-bind="value:replyNumber2"
                                onChange={(e) => this.setState({max_analise: e.target.value})}
                            />
                            <label for="text" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white">Доступно записей для поиска фактов EBM</label>
                            <input
                                type="number"
                                id="replyNumber3"
                                min="0"
                                step="1"
                                value={max_ddi}
                                data-bind="value:replyNumber3"
                                onChange={(e) => this.setState({max_ddi: e.target.value})}
                            />
                            <label for="text" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white">Dublicate Password</label>
                            <Select
                                className="basic-single"
                                classNamePrefix="select"
                                value={allow_type}
                                isSearchable
                                placeholder="Выберите класс"
                                name="topic"
                                options={allow_types}
                                getOptionLabel={(option) => option.label}
                                getOptionValue={(option) => option.value}
                                onChange={(x) => this.setState({allow_type: x})}
                            />
                            <br />
                            {message ?
                                messageStatus > 299 ?
                                  <p class="pb-2 mb-3 border-bottom" style={{ color: 'red' }}>{message}.</p>
                                  : messageStatus === 201 ?
                                    <p class="pb-2 mb-3 border-bottom" style={{ color: 'green' }}>{message}.</p>
                                  :null
                            : null}
                            <button type="button"
                                className="flex text-white bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:ring-blue-300 font-medium rounded-lg text-sm px-5 py-2.5 mr-2 dark:bg-blue-600 dark:hover:bg-blue-700 focus:outline-none dark:focus:ring-blue-800"
                                onClick={() => this.updateUserClick()}
                            >Обновить</button>
                        </div>
                    </div>

                </div>
            </div>
          </div>
        </>
      )
    }
  }
}